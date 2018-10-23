/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.core.LinearSolvers.SparseMatrix;
import flowpro.core.Parameters;
import java.io.IOException;

/**
 *
 * @author obublik
 */
abstract public class Preconditioner {
    
    
    public static Preconditioner factory(Parameters par) throws IOException {
        Preconditioner M = null;
        try {
            switch (par.preconditioner) {
                case "none":
                    M = new None(par);
                    break;
                case "jacobi":
                    M = new Jacobi(par);
                    break;
                case "blockjacobi":
                    M = new blockJacobi(par);
                    break;
                case "blockjacobiinversion":
                    M = new BlockJacobiInversion(par);
                    break;
                case "SSOR":
                    M = new SSOR(par);
                    break;
                case "ilu0":
                    M = new ILU0(par);
                    break;

                default:
                    throw new IOException("unknown preconditioner " + par.preconditioner);
            }
        } catch (Exception e) {
            System.out.println("Solver not set!");
        }

        return M;
    }
    
    abstract public void setMatrix(SparseMatrix A);
    
    abstract public void factor();
    
    abstract public void apply(double[] x, double[] b);
}