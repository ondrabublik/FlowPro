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
                    System.out.println("preconditioner: none");
                    break;
                case "jacobi":
                    M = new Jacobi(par);
                    System.out.println("preconditioner: jacobi");
                    break;
                case "blockjacobi":
                    M = new blockJacobi(par);
                    System.out.println("preconditioner: block jacobi");
                    break;
                case "blockjacobiinversion":
                    M = new BlockJacobiInversion(par);
                    System.out.println("preconditioner: inverted block jacobi");
                    break;
                case "ssor":
                    M = new SSOR(par);
                    System.out.println("preconditioner: SSOR");
                    break;
                case "ilu0":
                    M = new ILU0(par);
                    System.out.println("preconditioner: ILU(0)");
                    break;
                case "blockilu0":
                    M = new BlockILU0(par);
                    System.out.println("preconditioner: block ILU(0)");
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