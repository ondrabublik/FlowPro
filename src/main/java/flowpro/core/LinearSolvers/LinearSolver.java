package flowpro.core.LinearSolvers;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;
import flowpro.core.Parameters;
import flowpro.core.element.Element;
import java.io.IOException;
import java.util.Arrays;

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

    public enum LinearSolverType {
        gmres, umfpack, eigen, matlab;

        public static void help() {
            System.out.println("********************************");
            System.out.println("HELP for parameter linearSolver");
            System.out.println("list of possible values:");
            System.out.println(Arrays.asList(LinearSolverType.values()));
            System.out.println("********************************");
        }
    }
    
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
            LinearSolverType linearSolverType = LinearSolverType.valueOf(par.linearSolver.toLowerCase());
            switch (linearSolverType) {
                case gmres:
                    solver = new Gmres(A, M, 30, 5, par.iterativeSolverTol, par.nThreads);
                    System.out.println("linear solver: GMRES");
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

                case umfpack:
                    solver = new UmfpackExternSolver(A, par);
                    System.out.println("linear solver: UMFPACK");
                    break;
                    
                case eigen:
                    solver = new EigenExternSolver(A, par);
                    System.out.println("linear solver: EIGEN");
                    break;

                case matlab:
                    solver = new Matlab(A, par);
                    System.out.println("linear solver: MATLAB");
                    break;
            }
        } catch (IllegalArgumentException ex) {
            LinearSolverType.help();
            throw new IOException("unknown linear solver " + par.linearSolver.toLowerCase());
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
