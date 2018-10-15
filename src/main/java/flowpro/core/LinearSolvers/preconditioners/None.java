/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.core.LinearSolvers.SparseMatrixCRS;
import flowpro.core.Parameters;

/**
 *
 * @author obublik
 */
class None extends Preconditioner {

    None(Parameters par) {

    }

    @Override
    public void setMatrix(SparseMatrixCRS A) {
    }

    @Override
    public void factor() {
    }

    @Override
    public void apply(double[] x, double[] b) {
        System.arraycopy(b, 0, x, 0, b.length);
    }
}
