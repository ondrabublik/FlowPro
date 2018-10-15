/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.core.LinearSolvers.SparseMatrixCRS;
import flowpro.core.Parameters;
import java.util.Arrays;

/**
 * SSOR preconditioner. Uses symmetrical sucessive overrelaxation as a
 * preconditioner. Meant for symmetrical, positive definite matrices. For best
 * performance, omega must be carefully chosen (between 0 and 2).
 */
public class SSOR extends Preconditioner {

    double omegaF; // Overrelaxation parameter for the forward sweep
    double omegaR; // Overrelaxation parameter for the backwards sweep

    SparseMatrixCRS A;
    int n;
    private double[] xx; // Temporary vector for holding the half-step state

    /**
     * True if the reverse (backward) sweep is to be done. Without this, the
     * method is SOR instead of SSOR
     */
    private final boolean reverse;
    
    SSOR(Parameters par) {
        reverse = true;
        omegaF = 1;
        omegaR = 1;
    }

    @Override
    public void setMatrix(SparseMatrixCRS A) {
        this.A = A;
        n = A.getDofs();
        xx = new double[n];
    }

    @Override
    public void factor() {   
    }
    
    @Override
    public void apply(double[] xd, double[] bd) {
        int[] rowptr = A.getRowIndexes();
        int[] colind = A.getColumnIndexes();
        int[] diagind = A.getDiagonalIndexes();
        double[] data = A.getData();
        Arrays.fill(xd,0);
        System.arraycopy(xd, 0, xx, 0, n);

        // Forward sweep (xd oldest, xx halfiterate)
        for (int i = 0; i < n; ++i) {

            double sigma = 0;
            for (int j = rowptr[i]; j < diagind[i]; ++j)
                sigma += data[j] * xx[colind[j]];

            for (int j = diagind[i] + 1; j < rowptr[i + 1]; ++j)
                sigma += data[j] * xd[colind[j]];

            sigma = (bd[i] - sigma) / data[diagind[i]];

            xx[i] = xd[i] + omegaF * (sigma - xd[i]);
        }

        // Stop here if the reverse sweep was not requested
        if (!reverse) {
            System.arraycopy(xx, 0, xd, 0, n);
            return ;
        }

        // Backward sweep (xx oldest, xd halfiterate)
        for (int i = n - 1; i >= 0; --i) {

            double sigma = 0;
            for (int j = rowptr[i]; j < diagind[i]; ++j)
                sigma += data[j] * xx[colind[j]];

            for (int j = diagind[i] + 1; j < rowptr[i + 1]; ++j)
                sigma += data[j] * xd[colind[j]];

            sigma = (bd[i] - sigma) / data[diagind[i]];

            xd[i] = xx[i] + omegaR * (sigma - xx[i]);
        }
    }
}
