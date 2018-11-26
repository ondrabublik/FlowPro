package flowpro.core.DistributedLinearSolver;

import flowpro.core.DistributedLinearSolver.preconditioner.ParallelBlockJacobi;
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
    double[] x;
    double[][] V;
    double[] w, r, aux, b;
    int[] Imetis;

    int dofs, nThreads;
    int nLoad, nSave;
    int index;

    ParallelSparseMatrix A;
    ParallelBlockJacobi M;

    //public ParallelGmresSlave(SparseMatrix A, Preconditioner M, int m, int iterationLimit, double tol, int nThreads) {
    public ParallelGmresSlave(Element[] elems, int m, Parameters par, double[] x) throws IOException {

        this.elems = elems;
        this.x = x;
        nThreads = par.nThreads;

        dofs = 0;
        for (int i = 0; i < elems.length; i++) {
            Element elem = elems[i];
            int n = elem.nBasis * elem.getNEqs();
            dofs += n;
        }
        w = new double[dofs];
        r = new double[dofs];
        aux = new double[dofs];
        b = new double[dofs];
        V = new double[m + 1][dofs];

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

        // define sparse matrix
        A = new ParallelSparseMatrix(elems);
        Imetis = A.getImetis();
        M = new ParallelBlockJacobi(par);
        M.set(elems);
    }

    public MPIMessage doWork(MPIMessage msg) {
        MPIMessage msgOut = null;
        switch (msg.subTag) {
            case ParallelTags.UPDATE_MATRIXES_AND_VECTORS:
                // update A
                A.updateData();
                M.factor();

                // update b
                for (Element elem : metisElems) {
                    int[] ind = elem.gi_U;
                    for (int i = 0; i < ind.length; i++) {
                        b[ind[i]] = elem.RHS_loc[i];
                    }
                }

                msgOut = new MPIMessage(Tag.GMRES2MASTER, norm('b'));
                break;

            case ParallelTags.SMULT: // substract multiply
                A.SubstrMult(aux, b, x, nThreads);
                M.apply(r,aux);

                msgOut = new MPIMessage(Tag.GMRES2MASTER, norm('r'));
                break;

            case ParallelTags.MULT: // multiply
                index = (int) msg.getData();
                A.Mult(aux, V[index], nThreads);
                M.apply(w,aux);

                double[] innerProd = scalarProducts(index + 1);

                msgOut = new MPIMessage(Tag.GMRES2MASTER, innerProd);
                break;

            case ParallelTags.FILLV:
                double[] data = (double[]) msg.getData();
                int p = (int) data[0];
                double norm = data[1];
                double[] Vaux = V[p];
                if (p == 0) {
                    for (int i : Imetis) {
                        Vaux[i] = r[i] / norm;
                    }
                } else {
                    for (int i : Imetis) {
                        Vaux[i] = w[i] / norm;
                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;

            case ParallelTags.UPDATE_W:
                double[] H = (double[]) msg.getData();
                for (int k = 0; k < H.length; k++) {
                    Vaux = V[k];
                    for (int i : Imetis) {
                        w[i] -= Vaux[i] * H[k];
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
                            int[] ind = elem.gi_U;
                            double[] yElem = new double[ind.length];
                            for (int i = 0; i < ind.length; i++) {
                                yElem[i] = V[index][ind[i]];
                            }
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
                    int[] ind = elems[j].gi_U;
                    for (int i = 0; i < ind.length; i++) {
                        V[index][ind[i]] = dataRec.y[i];
                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;

            case ParallelTags.UPDATE_SOLUTION:
                Update up = (Update) msg.getData();
                double[] y = up.y;
                int n = up.i;
                for (int k = 0; k <= n; k++) {
                    Vaux = V[k];
                    for (int i : Imetis) {
                        x[i] += Vaux[i] * y[k];
                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;
        }
        return msgOut;
    }

    double[] scalarProducts(int i) {
        double[] sProd = new double[i];
        for (int k = 0; k < i; k++) {
            double[] Vaux = V[k];
            for (int j : Imetis) {
                sProd[k] += Vaux[j] * w[j];
            }
        }

        return sProd;
    }

    double norm(char s) {
        double n = 0;
        switch (s) {
            case 'b':
                for (int i : Imetis) {
                    n += b[i] * b[i];
                }
                break;

            case 'r':
                for (int i : Imetis) {
                    n += r[i] * r[i];
                }
                break;

            case 'w':
                for (int i : Imetis) {
                    n += w[i] * w[i];
                }
                break;
        }
        return n;
    }
}
