package flowpro.core.LinearSolvers;

import flowpro.api.Mat;
import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import flowpro.core.parallel.LiteElement;
import flowpro.core.parallel.Tag;
import java.io.IOException;
import litempi.MPIMessage;

/**
 *
 * @author obublik
 */
public class ParallelGmresSlave {

    Element[] elems;
    Element[] metisElems;
    GMRESdata[] gmDat;
    double[] x;

    int dofs, nThreads;
    int nLoad, nSave;
    int index;

    //public ParallelGmresSlave(SparseMatrix A, Preconditioner M, int m, int iterationLimit, double tol, int nThreads) {
    public ParallelGmresSlave(Element[] elems, int m, Parameters par, double[] x) throws IOException {

        this.elems = elems;
        this.x = x;

        gmDat = new GMRESdata[elems.length];
        dofs = 0;
        for (int i = 0; i < elems.length; i++) {
            Element elem = elems[i];
            int n = elem.nBasis * elem.getNEqs();
            dofs += n;
            if (elem.insideComputeDomain) {
                gmDat[i] = new GMRESdata(m, n);
            }
        }

        // number of saved elements
        int nMetis = 0;
        nSave = 0;
        for (Element elem : elems) {
            if (elem.gmresSave) {
                nSave++;
            }
            if (elem.insideMetisDomain) {
                nMetis++;
            }
        }

        // define elements inside metis domain
        int s = 0;
        metisElems = new Element[nMetis];
        for (Element elem : elems) {
            if (elem.insideMetisDomain) {
                metisElems[s] = elem;
                s++;
            }
        }
    }

    public MPIMessage doWork(MPIMessage msg) {
        MPIMessage msgOut = null;
        switch (msg.subTag) {
            case ParallelTags.UPDATE_PRECONDITIONER:
                // block diagonal preconditioner
                for (Element elem : metisElems) {
                    elem.Aux = Mat.invert(elem.ADiag);
                }

                msgOut = new MPIMessage(Tag.GMRES2MASTER, norm('b'));
                break;

            case ParallelTags.SMULT: // substract multiply
                for (Element elem : metisElems) {
                    int e = elem.index;
                    int[] ind = elem.gi_U;
                    double[][] A = elem.ADiag;
                    double[] b = elem.RHS_loc;
                    double[] aux = gmDat[e].aux;
                    for (int i = 0; i < ind.length; i++) {
                        aux[i] = b[i];
                        for (int j = 0; j < ind.length; j++) {
                            aux[i] -= A[j][i] * x[ind[j]];
                        }
                    }
                    for (int k = 0; k < elem.nFaces; k++) {
                        if (elem.TT[k] > -1) {
                            int[] indR = elems[elem.TT[k]].gi_U;
                            A = elem.ANeighs[k].MR;
                            for (int i = 0; i < ind.length; i++) {
                                for (int j = 0; j < indR.length; j++) {
                                    aux[i] -= A[j][i] * x[indR[j]];
                                }
                            }
                        }
                    }
                }
                applyPreconditioner('r');

                msgOut = new MPIMessage(Tag.GMRES2MASTER, norm('r'));
                break;

            case ParallelTags.MULT: // multiply
                index = (int) msg.getData();
                for (Element elem : metisElems) {
                    int e = elem.index;
                    double[][] A = elem.ADiag;
                    double[] V = gmDat[e].V[index];
                    double[] aux = gmDat[e].aux;
                    for (int i = 0; i < aux.length; i++) {
                        aux[i] = 0;
                        for (int j = 0; j < aux.length; j++) {
                            aux[i] += A[j][i] * V[j];
                        }
                    }
                    for (int k = 0; k < elem.nFaces; k++) {
                        if (elem.TT[k] > -1) {
                            A = elem.ANeighs[k].MR;
                            V = gmDat[elem.TT[k]].V[index];
                            for (int i = 0; i < A[0].length; i++) {
                                for (int j = 0; j < A.length; j++) {
                                    aux[i] += A[j][i] * V[j];
                                }
                            }
                        }
                    }

                }

                applyPreconditioner('w');
                
                double[] innerProd = scalarProducts(index + 1);
                
                msgOut = new MPIMessage(Tag.GMRES2MASTER, innerProd);
                break;

            case ParallelTags.FILLV:
                double[] data = (double[]) msg.getData();
                int p = (int) data[0];
                double norm = data[1];
                if (p == 0) {
                    for (Element elem : metisElems) {
                        int e = elem.index;
                        double[] V = gmDat[e].V[p];
                        double[] r = gmDat[e].r;
                        for (int j = 0; j < V.length; j++) {
                            V[j] = r[j] / norm;
                        }

                    }
                } else {
                    for (Element elem : metisElems) {
                        int e = elem.index;
                        double[] V = gmDat[e].V[p];
                        double[] w = gmDat[e].w;
                        for (int j = 0; j < V.length; j++) {
                            V[j] = w[j] / norm;
                        }

                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;

            case ParallelTags.UPDATE_W:
                double[] H = (double[]) msg.getData();
                for (Element elem : metisElems) {
                    int e = elem.index;
                    double[] w = gmDat[e].w;
                    for (int k = 0; k < H.length; k++) {
                        double[] V = gmDat[e].V[k];
                        for (int j = 0; j < V.length; j++) {
                            w[j] -= V[j] * H[k];
                        }
                    }

                }
                
                msgOut = new MPIMessage(Tag.GMRES2MASTER, norm('w'));
                break;

            case ParallelTags.NORM:
                msgOut = new MPIMessage(Tag.GMRES2MASTER, norm((char) msg.getData()));
                break;

            case ParallelTags.SAVE_DATA:
                index = (int) msg.getData();
                LiteElement[] dataSend = new LiteElement[nSave];
                int s = 0;
                if (index == -1) {
                    for (Element elem : elems) {
                        if (elem.gmresSave) {
                            int[] gind = elem.gi_U;
                            double[] yElem = new double[gind.length];
                            for (int i = 0; i < gind.length; i++) {
                                yElem[i] = x[gind[i]];
                            }
                            dataSend[s] = new LiteElement(elem.index, yElem);
                            s++;
                        }
                    }
                } else {
                    for (int e = 0; e < elems.length; e++) {
                        Element elem = elems[e];
                        if (elem.gmresSave) {
                            double[] V = gmDat[e].V[index];
                            double[] yElem = new double[V.length];
                            System.arraycopy(V, 0, yElem, 0, V.length);
                            dataSend[s] = new LiteElement(elem.index, yElem);
                            s++;
                        }
                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER, dataSend);
                break;

            case ParallelTags.LOAD_X:
                LiteElement[] dataReceive = (LiteElement[]) msg.getData();
                for (LiteElement dataRec : dataReceive) {
                    int j = dataRec.index;
                    int[] gind = elems[j].gi_U;
                    for (int i = 0; i < gind.length; i++) {
                        x[gind[i]] = dataRec.y[i];
                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;

            case ParallelTags.LOAD_W:
                dataReceive = (LiteElement[]) msg.getData();
                for (LiteElement dataRec : dataReceive) {
                    int j = dataRec.index;
                    double[] V = gmDat[j].V[index];
                    System.arraycopy(dataRec.y, 0, V, 0, V.length);
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;

            case ParallelTags.UPDATE_SOLUTION:
                Update up = (Update) msg.getData();
                double[] y = up.y;
                int n = up.i;
                for (Element elem : metisElems) {
                    int e = elem.index;
                    int[] ind = elem.gi_U;
                    for (int k = 0; k <= n; k++) {
                        double[] V = gmDat[e].V[k];
                        for (int j = 0; j < V.length; j++) {
                            x[ind[j]] += V[j] * y[k];
                        }
                    }

                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;
        }
        return msgOut;
    }

    void applyPreconditioner(char s) {
        switch (s) {
            case 'r':
                for (Element elem : metisElems) {
                    int e = elem.index;
                    double[][] A = elem.Aux;
                    double[] r = gmDat[e].r;
                    double[] aux = gmDat[e].aux;
                    for (int i = 0; i < r.length; i++) {
                        r[i] = 0;
                        for (int j = 0; j < r.length; j++) {
                            r[i] += A[j][i] * aux[j];
                        }
                    }
                }
                break;

            case 'w':
                for (Element elem : metisElems) {
                    int e = elem.index;
                    double[][] A = elem.Aux;
                    double[] w = gmDat[e].w;
                    double[] aux = gmDat[e].aux;
                    for (int i = 0; i < w.length; i++) {
                        w[i] = 0;
                        for (int j = 0; j < w.length; j++) {
                            w[i] += A[j][i] * aux[j];
                        }
                    }
                }
                break;
        }
    }

    double[] scalarProducts(int i) {
        double[] sProd = new double[i];
        for (Element elem : metisElems) {
            int e = elem.index;
            double[] w = gmDat[e].w;
            for (int k = 0; k < i; k++) {
                double[] V = gmDat[e].V[k];
                for (int j = 0; j < V.length; j++) {
                    sProd[k] += V[j] * w[j];
                }
            }
        }
        return sProd;
    }

    double norm(char s) {
        double n = 0;
        switch (s) {
            case 'b':
                for (Element elem : metisElems) {
                    double[] b = elem.RHS_loc;
                    for (int i = 0; i < b.length; i++) {
                        n += b[i] * b[i];
                    }
                }
                break;

            case 'r':
                for (Element elem : metisElems) {
                    int e = elem.index;
                    double[] r = gmDat[e].r;
                    for (int i = 0; i < r.length; i++) {
                        n += r[i] * r[i];
                    }
                }
                break;

            case 'w':
                for (Element elem : metisElems) {
                    int e = elem.index;
                    double[] w = gmDat[e].w;
                    for (int j = 0; j < w.length; j++) {
                        n += w[j] * w[j];
                    }
                }
                break;
        }
        return n;
    }
}
