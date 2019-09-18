/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.element.Element;
import flowpro.core.element.ImplicitBDFElement;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class SparseMatrixCRS {

    int dofs, nnz;
    public int[] IA, JA, diagind;
    public double[] HA;
    double[][][][] blocks;
    Element[] elems;

    SparseMatrixCRS(Element[] elems) {
        this.elems = elems;
        scanMatrixStructure();
    }

    private void scanMatrixStructure() {
        // alocate sparse matrix in crs
        dofs = computeDofs();
        nnz = computeNNZ();
        int[] IAaux = new int[nnz];
        JA = new int[nnz];
        HA = new double[nnz];

        IAaux[0] = 0;
        int si = 1;
        int sj = 0;
        blocks = new double[elems.length][][][];
        for (int i = 0; i < elems.length; i++) {
            Element elem = elems[i];
            int nBlockRow = 1; // number of blocks in a row
            for (int k = 0; k < elem.TT.length; k++) {
                if (elem.TT[k] > -1) {
                    nBlockRow++;
                }
            }

            blocks[i] = new double[nBlockRow][][]; // alocation of row matrixes
            int[] indexes = new int[nBlockRow];
            blocks[i][0] = ((ImplicitBDFElement)elem.ti).ADiag;
            indexes[0] = elem.index;
            nBlockRow = 1;
            for (int k = 0; k < elem.TT.length; k++) {
                if (elem.TT[k] > -1) {
                    blocks[i][nBlockRow] = ((ImplicitBDFElement)elem.ti).ANeighs[k].A;
                    indexes[nBlockRow] = elems[elem.TT[k]].index;
                    nBlockRow++;
                }
            }

            // sorting of matrix at rows acording to element indexes (bublesort)
            double[][] aux;
            int pom;
            for (int p = 0; p < nBlockRow - 1; p++) {
                for (int q = 0; q < nBlockRow - 1; q++) {
                    if (indexes[q] > indexes[q + 1]) {
                        aux = blocks[i][q];
                        blocks[i][q] = blocks[i][q + 1];
                        blocks[i][q + 1] = aux;
                        pom = indexes[q];
                        indexes[q] = indexes[q + 1];
                        indexes[q + 1] = pom;
                    }
                }
            }

            // create sparse structure
            int n = elem.getNEqs() * elem.nBasis;
            for (int k = 0; k < n; k++) {
                for (int q = 0; q < nBlockRow; q++) {
                    int[] glob = elems[indexes[q]].gi_U;
                    for (int j = 0; j < glob.length; j++) { // row cycle
                        JA[sj] = glob[j];
                        sj++;
                    }
                }
                IAaux[si] = sj;
                si++;
            }
        }

        IA = new int[si];
        System.arraycopy(IAaux, 0, IA, 0, si);
        IAaux = null;

        // Find the indices to the diagonal entries for ILU(0)
        diagind = findDiagonalIndices(dofs, JA, IA);
    }

    public int computeNNZ() {
        int s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            s += n * n;
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    s += n * ne;
                }
            }
        }
        return s;
    }

    public int computeDofs() {
        int s = 0;
        for (Element elem : elems) {
            s += elem.getNEqs() * elem.nBasis;
        }
        return s;
    }

    void updateData() {
        // fill A
        int s = 0;
        for (int i = 0; i < blocks.length; i++) {
            for (int p = 0; p < blocks[0][0][0].length; p++) {
                for (int j = 0; j < blocks[i].length; j++) {
                    double[][] A = blocks[i][j];
                    for (int k = 0; k < A.length; k++) { // row cycle
                        HA[s] = A[k][p];
                        s++;
                    }
                }
            }
        }
    }

    private int[] findDiagonalIndices(int m, int[] colind, int[] rowptr) {
        int[] diagI = new int[m];
        for (int k = 0; k < m; ++k) {
            diagI[k] = Arrays.binarySearch(colind, rowptr[k], rowptr[k + 1], k);
            if (diagI[k] < 0) {
                throw new RuntimeException("Missing diagonal entry on row " + (k + 1));
            }
        }

        return diagI;
    }

    public void updateB(double[] b) {
        int s = 0;
        for (Element elem : elems) {
            double[] RHS_loc = ((ImplicitBDFElement)elem.ti).RHS_loc;
            int n = elem.getNEqs() * elem.nBasis;
            for (int i = 0; i < n; i++) {
                b[s] = RHS_loc[i];
                s++;
            }
        }
    }

    public void Mult(double[] y, double[] x) {
        for (int i = 0; i < dofs; i++) {
            y[i] = 0;
            for (int j = IA[i]; j < IA[i + 1]; j++) {
                y[i] += HA[j] * x[JA[j]];
            }
        }
    }

    public void Mult(double[] y, double[] x, int nThreads) {
        MultThread[] parallel = new MultThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new MultThread(v, nThreads, x, y);
            parallel[v].start();
        }
        try {
            for (int v = 0; v < nThreads; v++) {
                parallel[v].join();
            }
        } catch (java.lang.InterruptedException e) {
            System.out.println(e);
        }
    }

    class MultThread extends Thread {

        int n, nStart, nThreads;
        double[] x, y;

        MultThread(int nStart, int nThreads, double[] x, double[] y) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.x = x;
            this.y = y;
            n = x.length;
        }

        @Override
        public void run() {
            for (int i = nStart; i < n; i = i + nThreads) {
                y[i] = 0;
                for (int j = IA[i]; j < IA[i + 1]; j++) {
                    y[i] += HA[j] * x[JA[j]];
                }
            }
        }
    }

    public void SubstrMult(double[] y, double[] b, double[] x) {
        for (int i = 0; i < dofs; i++) {
            y[i] = b[i];
            for (int j = IA[i]; j < IA[i + 1]; j++) {
                y[i] -= HA[j] * x[JA[j]];
            }
        }
    }

    public void SubstrMult(double[] y, double[] b, double[] x, int nThreads) {
        SubstrMultThread[] parallel = new SubstrMultThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new SubstrMultThread(v, nThreads, x, b, y);
            parallel[v].start();
        }
        try {
            for (int v = 0; v < nThreads; v++) {
                parallel[v].join();
            }
        } catch (java.lang.InterruptedException e) {
            System.out.println(e);
        }
    }

    class SubstrMultThread extends Thread {

        int n, nStart, nThreads;
        double[] x, b, y;

        SubstrMultThread(int nStart, int nThreads, double[] x, double[] b, double[] y) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.x = x;
            this.b = b;
            this.y = y;
            n = x.length;
        }

        @Override
        public void run() {
            for (int i = nStart; i < n; i = i + nThreads) {
                y[i] = b[i];
                for (int j = IA[i]; j < IA[i + 1]; j++) {
                    y[i] -= HA[j] * x[JA[j]];
                }
            }
        }
    }

    public int[] getRowIndexes() {
        return IA;
    }

    public int[] getColumnIndexes() {
        return JA;
    }

    public double[] getData() {
        return HA;
    }

    public int[] getDiagonalIndexes() {
        return diagind;
    }

    public int getNNZ() {
        return nnz;
    }

    public int getDofs() {
        return dofs;
    }

    public Element[] getElems() {
        return elems;
    }
}
