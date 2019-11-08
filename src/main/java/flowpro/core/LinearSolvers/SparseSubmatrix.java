/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.element.Element;
import flowpro.core.element.Implicit;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class SparseSubmatrix {

    int dofs, nnz;
    public int[] Icoo, Icrs, Iccs, Jcoo, Jcrs, Jccs, diagind, indexMap, globMap, globMapInv;
    public double[] H;
    Element[] elems;
    boolean[] inside;

    public SparseSubmatrix(Element[] elems, boolean[] inside) {
        this.elems = elems;
        this.inside = inside;
        defineGlobalIndexMap();
        scanMatrixStructure();
    }

    private void defineGlobalIndexMap() {
        int s = 0;
        int si = 0;
        for (int i = 0; i < elems.length; i++) {
            int n = elems[i].getNEqs() * elems[i].nBasis;
            s += n;
            if (inside[i]){
                si += n;
            }
        }
        globMap = new int[s];
        globMapInv = new int[si];
        
        s = 0;
        si = 0;
        for (int i = 0; i < elems.length; i++) {
            if (inside[i]){
                int[] gi = elems[i].gIndex;
                for(int j = 0; j < gi.length; j++){
                    globMap[s] = si;
                    globMapInv[si] = gi[j];
                    s++;
                    si++;
                }
            } else {
                s += elems[i].getNEqs() * elems[i].nBasis;
            }
        }
    }

    private void scanMatrixStructure() {
        dofs = computeDofs();
        nnz = computeNNZ();
        Icoo = new int[nnz];
        Jcoo = new int[nnz];
        H = new double[nnz];
        indexMap = new int[nnz];

        int s = 0;
        for (int r = 0; r < elems.length; r++) {
            if (inside[r]) {
                Element elem = elems[r];
                int n = elem.getNEqs() * elem.nBasis;
                int[] glob = elem.gIndex;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        Icoo[s] = globMap[glob[i]];
                        Jcoo[s] = globMap[glob[j]];
                        indexMap[s] = s;
                        s++;
                    }
                }
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TT[k] > -1 && inside[elems[elem.TT[k]].index]) {
                        int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                        int[] globe = elems[elem.TT[k]].gIndex;
                        for (int i = 0; i < n; i++) {
                            for (int j = 0; j < ne; j++) {
                                Icoo[s] = globMap[glob[i]];
                                Jcoo[s] = globMap[globe[j]];
                                indexMap[s] = s;
                                s++;
                            }
                        }
                    }
                }
            }
        }
    }

    public void buildCRSformat() {
        quickSort(Icoo, Jcoo, indexMap, 0, nnz - 1);
        indexMap = invertMap(indexMap);
        Icrs = new int[nnz];
        Jcrs = new int[nnz];
        System.arraycopy(Icoo, 0, Icrs, 0, nnz);
        System.arraycopy(Jcoo, 0, Jcrs, 0, nnz);
        Icrs = compress(Icrs);
        diagind = findDiagonalIndices(dofs, Jcrs, Icrs);
    }

    public int computeNNZ() {
        int s = 0;
        for (int i = 0; i < elems.length; i++) {
            if (inside[i]) {
                Element elem = elems[i];
                int n = elem.getNEqs() * elem.nBasis;
                s += n * n;
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TT[k] > -1 && inside[elems[elem.TT[k]].index]) {
                        int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                        s += n * ne;
                    }
                }
            }
        }
        return s;
    }

    public int computeDofs() {
        int s = 0;
        for (int i = 0; i < elems.length; i++) {
            if (inside[i]) {
                s += elems[i].getNEqs() * elems[i].nBasis;
            }
        }
        return s;
    }

    public void updateData() {
        int s = 0;
        for (int r = 0; r < elems.length; r++) {
            if (inside[r]) {
                Element elem = elems[r];
                int n = elem.getNEqs() * elem.nBasis;
                double[][] Ad = ((Implicit)elem.ti).ADiag;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        H[indexMap[s]] = Ad[j][i];
                        s++;
                    }
                }
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TT[k] > -1 && inside[elems[elem.TT[k]].index]) {
                        int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                        double[][] An = ((Implicit)elem.ti).ANeighs[k].A;
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
    }

    private int[] findDiagonalIndices(int m, int[] colind, int[] rowptr) {
        int[] diagI = new int[m];
        for (int k = 0; k < m; ++k) {
            if (rowptr[k] != rowptr[k + 1]) {
                diagI[k] = Arrays.binarySearch(colind, rowptr[k], rowptr[k + 1], k);
                if (diagI[k] < 0) {
                    throw new RuntimeException("Missing diagonal entry on row " + (k + 1));
                }
            }
        }
        return diagI;
    }

    public void Mult(double[] y, double[] x) {
        for (int i = 0; i < dofs; i++) {
            int ig = globMapInv[i];
            y[ig] = 0;
            for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                y[ig] += H[j] * x[Jcrs[j]];
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
                int ig = globMapInv[i];
                y[ig] = 0;
                for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                    y[ig] += H[j] * x[Jcrs[j]];
                }
            }
        }
    }

    public void SubstrMult(double[] y, double[] b, double[] x) {
        for (int i = 0; i < dofs; i++) {
            int ig = globMapInv[i];
            y[ig] = b[ig];
            for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                y[ig] -= H[j] * x[Jcrs[j]];
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
                int ig = globMapInv[i];
                y[ig] = b[ig];
                for (int j = Icrs[i]; j < Icrs[i + 1]; j++) {
                    y[ig] -= H[j] * x[Jcrs[j]];
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
        int[] Icomp = new int[dofs + 1];
        int s = 0;
        for (int i = 1; i < dofs; i++) {
            int sum = 0;
            while (s < nnz && I[s] == i - 1) {
                s++;
                sum++;

            }
            Icomp[i] = Icomp[i - 1] + sum;
        }
        Icomp[dofs] = nnz;

        return Icomp;
    }

    private void quickSort(int[] array, int[] array2, int[] index, int lowerIndex, int higherIndex) {
        int i = lowerIndex;
        int j = higherIndex;
        int pivotIndex = lowerIndex + (higherIndex - lowerIndex) / 2;
        int pivot = array[pivotIndex];
        int pivot2 = array2[pivotIndex];
        while (i <= j) {
            while (array[i] < pivot || (array[i] == pivot && array2[i] < pivot2)) {
                i++;
            }
            while (array[j] > pivot || (array[j] == pivot && array2[j] > pivot2)) {
                j--;
            }
            if (i <= j) {
                swap(array, array2, index, i, j);
                i++;
                j--;
            }
        }
        if (lowerIndex < j) {
            quickSort(array, array2, index, lowerIndex, j);
        }
        if (i < higherIndex) {
            quickSort(array, array2, index, i, higherIndex);
        }
    }

    private void quickSort(int[] array, int[] array2, int lowerIndex, int higherIndex) {
        int i = lowerIndex;
        int j = higherIndex;
        int pivot = array[lowerIndex + (higherIndex - lowerIndex) / 2];
        while (i <= j) {
            while (array[i] < pivot) {
                i++;
            }
            while (array[j] > pivot) {
                j--;
            }
            if (i <= j) {
                swap(array, array2, i, j);
                i++;
                j--;
            }
        }
        if (lowerIndex < j) {
            quickSort(array, array2, lowerIndex, j);
        }
        if (i < higherIndex) {
            quickSort(array, array2, i, higherIndex);
        }
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
        for (int i = 0; i < nnz; i++) {
            index2[i] = i;
        }
        quickSort(index, index2, 0, nnz - 1);
        return index2;
    }
}
