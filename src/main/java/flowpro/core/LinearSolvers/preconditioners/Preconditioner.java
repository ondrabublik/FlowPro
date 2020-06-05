/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.core.LinearSolvers.SparseMatrix;
import flowpro.core.Parameters;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
abstract public class Preconditioner {
    
    public enum PreconditionerType {
        none, jacobi, blockjacobi, blockjacobiinversion, ssor, ilu0, blockilu0;

        public static void help() {
            System.out.println("********************************");
            System.out.println("HELP for parameter preconditioner");
            System.out.println("list of possible values:");
            System.out.println(Arrays.asList(PreconditionerType.values()));
            System.out.println("********************************");
        }
    }
    
    public static Preconditioner factory(Parameters par) throws IOException {
        Preconditioner M = null;
        PreconditionerType preconditionerType = PreconditionerType.valueOf(par.preconditioner);
        try {
            switch (preconditionerType) {
                case none:
                    M = new None(par);
                    break;
                case jacobi:
                    M = new Jacobi(par);
                    break;
                case blockjacobi:
                    M = new blockJacobi(par);
                    break;
                case blockjacobiinversion:
                    M = new BlockJacobiInversion(par);
                    break;
                case ssor:
                    M = new SSOR(par);
                    break;
                case ilu0:
                    M = new ILU0(par);
                    break;
                case blockilu0:
                    M = new BlockILU0(par);
                    break;
            }
        } catch (IllegalArgumentException ex) {
            PreconditionerType.help();
            throw new IOException("unknown preconditioner " + par.preconditioner);
        }

        return M;
    }
    
    abstract public void setMatrix(SparseMatrix A);
    
    abstract public void factor();
    
    abstract public void apply(double[] x, double[] b);
}