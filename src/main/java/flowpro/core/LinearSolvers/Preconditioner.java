/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import java.util.Arrays;

/**
 *
 * @author obublik
 */
class Preconditioner {

    SparseMatrixCRS A;
    double[] ilu;
    int n;

    Preconditioner(SparseMatrixCRS A) {
        this.A = A;
        ilu = new double[A.getNNZ()];
        n = A.getDofs();
    }

    public void factor() {
        // Internal CRS matrix storage 
        int[] IA = A.getRowIndexes();
        int[] JA = A.getColumnIndexes();
        double[] HA = A.getData();
        int[] diagind = A.getDiagonalIndexes();
        System.arraycopy(HA, 0, ilu, 0, HA.length);

        // Go down along the main diagonal 
        for (int k = 1; k < n; ++k) {
            for (int i = IA[k]; i < diagind[k]; ++i) {

                // Get the current diagonal entry 
                int index = JA[i];
                double LUii = ilu[diagind[index]];

                if (LUii == 0) {
                    throw new RuntimeException("Zero pivot encountered on row "
                            + (i + 1) + " during ILU process");
                }

                // Elimination factor 
                double LUki = (ilu[i] /= LUii);

                // Traverse the sparse row i, reducing on row k 
                for (int j = diagind[index] + 1, l = IA[k] + 1; j < IA[index + 1]; ++j) {

                    while (l < IA[k + 1] && JA[l] < JA[j]) {
                        l++;
                    }

                    if (JA[l] == JA[j]) {
                        ilu[l] -= LUki * ilu[j];
                    }
                }
            }
        }
    }

    public void apply(double[] x, double[] b) {
        int[] IA = A.getRowIndexes();
        int[] JA = A.getColumnIndexes();
        int[] diagind = A.getDiagonalIndexes();
        double[] y = new double[x.length];
        Arrays.fill(x, 0);

        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = IA[i]; j < diagind[i]; j++) {
                sum = sum + ilu[j] * y[JA[j]];
            }
            y[i] = b[i] - sum;
        }

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = diagind[i] + 1; j < IA[i + 1]; j++) {
                sum = sum + ilu[j] * x[JA[j]];
            }
            x[i] = (y[i] - sum) / ilu[diagind[i]];
        }
    }

    public double[] getData() {
        return ilu;
    }
}
