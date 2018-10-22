/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.Mesh.Element;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class SparseMatrix {

    int dofs, nnz;
    public int[] Icoo, Icrs, Iccs, Jcoo, Jcrs, Jccs, diagind, indexMap;
    public double[] H;
    Element[] elems;

    SparseMatrix(Element[] elems) {
        this.elems = elems;
        scanMatrixStructure();
    }

    private void scanMatrixStructure() {
        dofs = computeDofs();
        nnz = computeNNZ();
        Icoo = new int[nnz];
        Jcoo = new int[nnz];
        H = new double[nnz];
        indexMap = new int[nnz];

        int s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            int[] glob = elem.gi_U;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Icoo[s] = glob[i];
                    Jcoo[s] = glob[j];
                    indexMap[s] = s;
                    s++;
                }
            }
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    int[] globe = elems[elem.TT[k]].gi_U;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < ne; j++) {
                            Icoo[s] = glob[i];
                            Jcoo[s] = globe[j];
                            indexMap[s] = s;
                            s++;
                        }
                    }
                }
            }
        }
    }

    public void buildCRSformat() {
        quicksort(Icoo, Jcoo, indexMap, 0, nnz);
        indexMap = invertMap(indexMap);
        Icrs = new int[nnz];
        Jcrs = new int[nnz];
        System.arraycopy(Icoo, 0, Icrs, 0, nnz);
        System.arraycopy(Jcoo, 0, Jcrs, 0, nnz);
        Icrs = compress(Icrs);
        diagind = findDiagonalIndices(dofs, Jcrs, Icrs);
    }

    public void buildCCSformat() {
        quicksort(Jcoo, Icoo, indexMap, 0, nnz);
        indexMap = invertMap(indexMap);
        Iccs = new int[nnz];
        Jccs = new int[nnz];
        System.arraycopy(Icoo, 0, Iccs, 0, nnz);
        System.arraycopy(Jcoo, 0, Jccs, 0, nnz);
        Jccs = compress(Jccs);
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
        int s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            double[][] Ad = elem.ADiag;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    H[indexMap[s]] = Ad[j][i];
                    s++;
                }
            }
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    double[][] An = elem.ANeighs[k].MR;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < ne; j++) {
                            H[indexMap[s]] = An[j][i];
                            s++;
                        }
                    }
                }
            }
        }
    }

    public void updateB(double[] b) {
        int s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            for (int i = 0; i < n; i++) {
                b[s] = elem.RHS_loc[i];
                s++;
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

    public void Mult(double[] y, double[] x) {
        for (int i = 0; i < dofs; i++) {
            y[i] = 0;
            for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                y[i] += H[j] * x[Jcrs[j]];
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
                for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                    y[i] += H[j] * x[Jcrs[j]];
                }
            }
        }
    }

    public void SubstrMult(double[] y, double[] b, double[] x) {
        for (int i = 0; i < dofs; i++) {
            y[i] = b[i];
            for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                y[i] -= H[j] * x[Jcrs[j]];
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
                for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                    y[i] -= H[j] * x[Jcrs[j]];
                }
            }
        }
    }

    public int[] getRowIndexesCRS() {
        return Icrs;
    }

    public int[] getColumnIndexesCRS() {
        return Jcrs;
    }

    public int[] getRowIndexesCCS() {
        return Iccs;
    }

    public int[] getColumnIndexesCCS() {
        return Jccs;
    }

    public int[] getRowIndexesCOO() {
        return Icoo;
    }

    public int[] getColumnIndexesCOO() {
        return Jcoo;
    }

    public double[] getData() {
        return H;
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

    public int[] compress(int[] I) {
        int[] J = new int[dofs + 1];
        int s = 0;
        for (int i = 0; i < nnz; i++) {
            if (I[i] != s) {
                s++;
                J[s] = i;
            }
        }
        J[dofs] = nnz;
        return J;
    }

    public static void quicksort(int[] array, int[] array2, int[] index, int left, int right) {
        if (left < right) {
            int boundary = left;
            for (int i = left + 1; i < right; i++) {
                if (array[i] < array[left] || (array[i] == array[left] && array2[i] < array2[left])) {
                    swap(array, array2, index, i, ++boundary);
                }
            }
            swap(array, array2, index, left, boundary);
            quicksort(array, array2, index, left, boundary);
            quicksort(array, array2, index, boundary + 1, right);
        }
    }

    public static void quicksort(int[] array, int[] array2, int left, int right) {
        if (left < right) {
            int boundary = left;
            for (int i = left + 1; i < right; i++) {
                if (array[i] < array[left]) {
                    swap(array, array2, i, ++boundary);
                }
            }
            swap(array, array2, left, boundary);
            quicksort(array, array2, left, boundary);
            quicksort(array, array2, boundary + 1, right);
        }
    }

    private static void swap(int[] array, int[] array2, int[] index, int left, int right) {
        int tmp = array[right];
        array[right] = array[left];
        array[left] = tmp;

        tmp = array2[right];
        array2[right] = array2[left];
        array2[left] = tmp;

        tmp = index[right];
        index[right] = index[left];
        index[left] = tmp;
    }

    private static void swap(int[] array, int[] array2, int left, int right) {
        int tmp = array[right];
        array[right] = array[left];
        array[left] = tmp;

        tmp = array2[right];
        array2[right] = array2[left];
        array2[left] = tmp;
    }

    public int[] invertMap(int[] index) {
        int[] index2 = new int[nnz];
        for(int i = 0; i < nnz; i++){
            index2[i] = i;
        }
        quicksort(index, index2, 0, nnz);
        return index2;
    }
}
