/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.DistributedLinearSolver.preconditioner;

import flowpro.core.LinearSolvers.SparseSubmatrix;
import flowpro.core.Parameters;
import flowpro.core.element.Element;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class ParallelIlu0 extends ParallelPreconditioner {

    Element[] elems;
    SparseSubmatrix A;
    int n;
    double[] ilu;

    ParallelIlu0(Parameters par) {
        System.out.println("Preconditioner: ilu0");
    }

    public void set(Element[] elems) {
        this.elems = elems;
        boolean[] inside = new boolean[elems.length];
        for (int i = 0; i < elems.length; i++) {
            inside[i] = elems[i].insideComputeDomain;
        }
        A = new SparseSubmatrix(elems, inside);
        A.buildCRSformat();
        ilu = new double[A.getNNZ()];
        n = A.getDofs();
    }

    @Override
    public void factor() {
        // update data
        A.updateData();
        
        // Internal CRS matrix storage 
        int[] IA = A.getRowIndexesCRS();
        int[] JA = A.getColumnIndexesCRS();
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

    @Override
    public void apply(double[] x, double[] b) {
        int[] IA = A.getRowIndexesCRS();
        int[] JA = A.getColumnIndexesCRS();
        int[] diagind = A.getDiagonalIndexes();
        int[] globMapInv = A.globMapInv;
        double[] y = new double[x.length];
        Arrays.fill(x, 0);

        for (int i = 0; i < n; i++) {
            int gii = globMapInv[i];
            double sum = 0;
            for (int j = IA[i]; j < diagind[i]; j++) {
                sum = sum + ilu[j] * y[globMapInv[JA[j]]];
            }
            y[gii] = b[gii] - sum;
        }

        for (int i = n - 1; i >= 0; i--) {
            int gii = globMapInv[i];
            double sum = 0;
            for (int j = diagind[i] + 1; j < IA[i + 1]; j++) {
                sum = sum + ilu[j] * x[globMapInv[JA[j]]];
            }
            x[gii] = (y[gii] - sum) / ilu[diagind[i]];
        }
    }
}
