/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package flowpro.core.DistributedLinearSolver.preconditioner;

import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import java.io.IOException;

/**
 *
 * @author obublik
 */
abstract public class ParallelPreconditioner {
    
    public static ParallelPreconditioner factory(Parameters par) throws IOException {
        ParallelPreconditioner M = null;
        try {
            switch (par.preconditioner) {
                case "blockjacobiinversion":
                    M = new ParallelBlockJacobiInversion(par);
                    break;
                    
                case "ilu0":
                    M = new ParallelIlu0(par);
                    break;

                default:
                    throw new IOException("unknown preconditioner " + par.preconditioner);
            }
        } catch (Exception e) {
            System.out.println("Solver not set!");
        }

        return M;
    }
    
    abstract public void set(Element[] elems);
    
    abstract public void factor();
    
    abstract public void apply(double[] x, double[] b);
}