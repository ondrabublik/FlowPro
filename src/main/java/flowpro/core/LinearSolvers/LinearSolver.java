package flowpro.core.LinearSolvers;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;
import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import java.io.IOException;

/**
 *
 * @author obublik
 */
abstract public class LinearSolver {

    public static SparseMatrix A;
    public static Preconditioner M;
    public static double[] b;

    public static LinearSolver solver;

    abstract public boolean solve(double[] x, double[] b);

    public static LinearSolver factory(Element[] elems, Parameters par) throws IOException {

        // build matrix structure
        A = new SparseMatrix(elems);
        
        // define preconditiner
        M = Preconditioner.factory(par);
        M.setMatrix(A);

        // alocate RHS
        b = new double[A.getDofs()];
        
        solver = null;
        
        try {
            switch (par.linearSolver.toLowerCase()) {
                case "gmres":
                    solver = new Gmres(A, M, 30, 5, par.iterativeSolverTol, par.nThreads);
                    break;
//                case "jacobi":
//                    solver = new Jacobi(elems, dofs, 500, par.iterativeSolverTol, par.nThreads);
//                    break;
//
//                case "bicgstab":
//                    solver = new BiCgStab(elems, dofs, 100, par.iterativeSolverTol, par.nThreads);
//                    break;
//
//                case "extern":
//                    solver = new ExternSolver(elems, dofs, par);
//                    break;

                case "umfpack":
                    solver = new UmfpackExternSolver(A, par);
                    break;
                    
                case "eigen":
                    solver = new EigenExternSolver(A, par);
                    break;

                case "matlab":
                    solver = new Matlab(A, par);
                    break;   
                    
                default:
                    throw new IOException("unknown solver " + par.linearSolver);
            }
        } catch (Exception e) {
            System.out.println("Solver not set!");
        }

        return solver;
    }

    public boolean solve(double[] x) {

        A.updateData(); // update data in matrix A
        A.updateB(b);   // update RHS
        M.factor();     // update preconditioner M
        
        boolean flag = solver.solve(x, b);  // solve
        return flag;
    }
}
