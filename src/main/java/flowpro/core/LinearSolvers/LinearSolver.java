package flowpro.core.LinearSolvers;

import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import java.io.IOException;

/**
 *
 * @author obublik
 */
abstract public class LinearSolver {

    public Element[] elems;

    abstract public boolean solve(double[] x);

    public static LinearSolver factory(Parameters par, Element[] elems, int dofs) throws IOException {
        LinearSolver solver = null;
        try {
            switch (par.linearSolver) {
                case "jacobi":
                    solver = new Jacobi(elems, dofs, 500, par.iterativeSolverTol, par.nThreads);
                    break;

                case "bicgstab":
                    solver = new BiCgStab(elems, dofs, 100, par.iterativeSolverTol, par.nThreads);
                    break;

                case "extern":
                    solver = new ExternSolver(elems, dofs, par);
                    break;
                    
                case "externcppsolver":
                    solver = new ExternCppSolver(elems, dofs, par);
                    break;
                
                case "MTJ":
                    //solver = new MTJsolver(elems, dofs, 500, par.iterativeSolverTol, par.nThreads);
                    break;
                    
                case "new":
                    solver = new NewLinSol(elems, dofs, par);
                    break;   
                    
                default:
                    solver = new Gmres(elems, dofs, 30, 5, par.iterativeSolverTol, par.nThreads);
                    break;
            }
        } catch (Exception e) {
            System.out.println("Solver not set!");
        }

        return solver;
    }

    double scalarProduct(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            s = s + a[i] * b[i];
        }
        return s;
    }

    double[] copy(double[] a) {
        double[] b = new double[a.length];
        System.arraycopy(a, 0, b, 0, a.length);

        return b;
    }
}
