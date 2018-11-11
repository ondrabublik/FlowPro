package flowpro.core.LinearSolvers;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;
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
    SparseMatrix A;
    SparseMatrixBlock Ablock;
    Preconditioner M;
    double[] b, x;

    int dofs, nThreads;
    double[] r, aux;

    //public ParallelGmresSlave(SparseMatrix A, Preconditioner M, int m, int iterationLimit, double tol, int nThreads) {
    public ParallelGmresSlave(Element[] elems, Parameters par, double[] x) throws IOException {

        this.elems = elems;
        this.x = x;

        // build matrix structure
        A = new SparseMatrix(elems);
        dofs = A.dofs;
        b = new double[dofs];
        A.updateB(b);

        // define preconditiner
        //M = Preconditioner.factory(par);
        //M.setMatrix(A);

        // define domain block
        int sIn = 0;
        for (int i = 0; i < elems.length; i++) {
            if (elems[i].insideMetisDomain) {
                sIn += elems[i].nBasis * elems[i].getNEqs();
            }
        }
        int[] Iblock = new int[sIn];
        int s = 0;
        for (int i = 0; i < elems.length; i++) {
            if (elems[i].insideMetisDomain) {
                for (int j = 0; j < elems[i].nBasis * elems[i].getNEqs(); j++) {
                    Iblock[s] = elems[i].gi_U[j];
                    s++;
                }
            }
        }
        int[] Jblock = new int[dofs];
        for (int i = 0; i < dofs; i++) {
            Jblock[i] = i;
        }

        Ablock = new SparseMatrixBlock(Iblock, Jblock, A);

        // alocate RHS
        aux = new double[Ablock.getNRows()];
    }

    public MPIMessage doWork(MPIMessage msg) {
        MPIMessage msgOut = null;
        switch (msg.subTag) {
            case 0:
                //Ablock.SubstrMult(aux, b, x);
                //M.apply(r, aux);
                Ablock.SubstrMult(r, b, x);
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;
            case 1:
                for (int j = 0; j < n; j++) {
                    x[j] += r[j];
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;
            case 2:
                double error = norm(aux);
                msgOut = new MPIMessage(Tag.GMRES2MASTER, error);
                break;
            case 3:
                int nSave = Ablock.nIBlock;
                LiteElement[] dataSend = new LiteElement[nSave];
                for (int i = 0; i < nSave; i++) {
                    int j = Ablock.Iblock[i];
                    int nBasis = elems[j].nBasis;
                    int nEqs = elems[j].getNEqs();
                    double[] yElem = new double[nEqs * nBasis];
                    for (int m = 0; m < nEqs; m++) {
                        for (int p = 0; p < nBasis; p++) {
                            int ind = nBasis * m + p;
                            yElem[ind] = x[elems[j].gi_U[ind]];
                        }
                    }
                    dataSend[i] = new LiteElement(j, yElem);
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER, dataSend);
                break;
            case 4:
                LiteElement[] dataReceive = (LiteElement[]) msg.getData();
                for (LiteElement dataReceive1 : dataReceive) {
                    int j = dataReceive1.index;
                    for (int m = 0; m < nEqs; m++) {
                        for (int p = 0; p < elems[j].nBasis; p++) {
                            int ind = elems[j].nBasis * m + p;
                            x[elems[j].gi_U[ind]] = dataReceive1.y[ind];
                        }
                    }
                }
                msgOut = new MPIMessage(Tag.GMRES2MASTER);
                break;
        }
        return msgOut;
    }

    double norm(double[] a) {
        double n = 0;
        for (int i = 0; i < a.length; i++) {
            n = n + a[i] * a[i];
        }
        n = Math.sqrt(n);

        return n;
    }
}
