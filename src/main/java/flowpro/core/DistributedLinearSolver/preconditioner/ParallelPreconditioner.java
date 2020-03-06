/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package flowpro.core.DistributedLinearSolver.preconditioner;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;
import flowpro.core.Parameters;
import flowpro.core.element.Element;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
abstract public class ParallelPreconditioner {
    
    public enum ParallelPreconditionerType {
        blockjacobiinversion, ilu0;

        public static void help() {
            System.out.println("********************************");
            System.out.println("HELP for parameter preconditioner");
            System.out.println("list of possible values:");
            System.out.println(Arrays.asList(ParallelPreconditionerType.values()));
            System.out.println("********************************");
        }
    }
    
    public static ParallelPreconditioner factory(Parameters par) throws IOException {
        ParallelPreconditioner M = null;
        ParallelPreconditionerType preconditionerType = ParallelPreconditionerType.valueOf(par.parallelPreconditioner);
        try {
            switch (preconditionerType) {
                case blockjacobiinversion:
                    M = new ParallelBlockJacobiInversion(par);
                    break;
                    
                case ilu0:
                    M = new ParallelIlu0(par);
                    break;
            }
        } catch (IllegalArgumentException ex) {
            ParallelPreconditioner.ParallelPreconditionerType.help();
            throw new IOException("unknown preconditioner " + par.parallelPreconditioner);
        }

        return M;
    }
    
    abstract public void set(Element[] elems);
    
    abstract public void factor();
    
    abstract public void apply(double[] x, double[] b);
}