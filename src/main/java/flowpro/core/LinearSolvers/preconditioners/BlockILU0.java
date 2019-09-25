/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.api.Mat;
import flowpro.core.LinearSolvers.SparseMatrix;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import flowpro.core.element.Implicit;

/**
 *
 * @author obublik
 */
public class BlockILU0 extends Preconditioner {

    SparseMatrix A;
    int nThreads;
    int[][] indexes;
    double[][][] diagonalInverse;
    double[][][][] ILU;
    double[][] y;
    Element[] elems;
    int nElem;

    BlockILU0(Parameters par) {
        nThreads = par.nThreads;
    }

    @Override
    public void setMatrix(SparseMatrix A) {
        elems = A.getElems();
        nElem = elems.length;
        diagonalInverse = new double[nElem][][];
        ILU = new double[nElem][][][];
        indexes = new int[nElem][];
        for (int i = 0; i < nElem; i++) {
            Element elem = elems[i];
            int nNeigh = 1;
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    nNeigh++;
                }
            }
            ILU[i] = new double[nNeigh][][];
            indexes[i] = new int[nNeigh];
            // diagonal block
            int s = 0;
            int ne = elem.nBasis * elem.getNEqs();
            ILU[i][s] = new double[ne][ne];
            double[][] AuxILU = ILU[i][s];
            double[][] AuxA = ((Implicit)elem.ti).ADiag;
            for (int n = 0; n < ne; n++) {
                System.arraycopy(AuxA[n], 0, AuxILU[n], 0, ne);
            }
            indexes[i][s] = elem.index;
            s++;
            // out of diagonal block
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int me = elems[elem.TT[k]].nBasis * elems[elem.TT[k]].getNEqs();
                    ILU[i][s] = new double[ne][me];
                    AuxILU = ILU[i][s];
                    AuxA = ((Implicit)elem.ti).ANeighs[k].A;
                    for (int n = 0; n < ne; n++) {
                        System.arraycopy(AuxA[n], 0, AuxILU[n], 0, me);
                    }
                    indexes[i][s] = elems[elem.TT[k]].index;
                    s++;
                }
            }
        }

        // sorting
        for (int i = 0; i < nElem; i++) {
            // bubble sort
            for (int j = 0; j < indexes[i].length - 1; j++) {
                for (int k = 0; k < indexes[i].length - 1; k++) {
                    if (indexes[i][k] > indexes[i][k + 1]) {
                        int pom = indexes[i][k];
                        indexes[i][k] = indexes[i][k + 1];
                        indexes[i][k + 1] = pom;

                        double[][] Apom = ILU[i][k];
                        ILU[i][k] = ILU[i][k + 1];
                        ILU[i][k + 1] = Apom;
                    }
                }
            }
        }

        y = new double[nElem][];
        for (int i = 0; i < nElem; i++) {
            y[i] = new double[elems[i].nBasis * elems[i].getNEqs()];
        }
    }

    @Override
    public void factor() {
        for (int i = 1; i < nElem; i++) {
            int[] row = indexes[i];
            for (int k = 0; k < row.length; k++) {
                double[][] iAkk = Mat.invert(ILU[row[k]][find(row[k], row[k])]);
                if (row[k] < i - 1) {
                    double[][] Aik = ILU[i][k];
                    Aik = Mat.times(Aik, iAkk);
                    for (int j = 0; j < row.length; j++) {
                        int kj = find(row[k], row[j]);
                        if (row[j] > row[k] && kj != -1) {
                            double[][] Aij = ILU[i][j];
                            double[][] Akj = ILU[row[k]][kj];
                            Aij = Mat.minus(Aij, Mat.times(Aik, Akj));
                        }
                    }
                }
            }
        }
    }

    @Override
    public void apply(double[] x, double[] b) {
        for (int i = 0; i < nElem; i++) {
            int[] row = indexes[i];
            double[] yp = y[i];
            int[] gi = elems[i].gIndex;
            for(int j = 0; j < gi.length; j++){
                yp[j] = b[gi[j]];
            }
            for (int k = 0; k < row.length; k++) {
                if (row[k] < i) {
                    double[][] L = ILU[i][k];
                    double[] yL = y[row[k]];
                    for (int s = 0; s < L.length; s++) {
                        for (int t = 0; t < L[0].length; t++) {
                            yp[s] -= L[s][t] * yL[t];
                        }
                    }
                }
            }
        }

        for (int i = nElem - 1; i >= 0; i--) {
            int[] row = indexes[i];
            double[] sum = new double[y[i].length];
            for (int k = 0; k < row.length; k++) {
                if (row[k] > i) {
                    int[] gi = elems[row[k]].gIndex;
                    double[][] U = ILU[i][k];
                    for (int s = 0; s < U.length; s++) {
                        for (int t = 0; t < U[0].length; t++) {
                            sum[s] += U[s][t] * x[gi[t]];
                        }
                    }
                }
            }
            
            double[][] iUii = Mat.invert(ILU[i][find(i, i)]);
            double[] yp = y[i];
            for(int j = 0; j < yp.length; j++){
                sum[j] = yp[j] - sum[j];
            }
            sum = Mat.times(iUii, sum);
            int[] gi = elems[i].gIndex;
            for(int j = 0; j < gi.length; j++){
                x[gi[j]] = sum[j];
            }
        }

//        for (int i = 0; i < b.length; i++) {
//            x[i] = b[i];
//        }
    }

    private int find(int row, int column) {
        int[] ind = indexes[row];
        for (int i = 0; i < ind.length; i++) {
            if (ind[i] == column) {
                return i;
            }
        }
        return -1;
    }
}
