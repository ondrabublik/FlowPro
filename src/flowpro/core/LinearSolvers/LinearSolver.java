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

    void ComputeResiduum(double[] x, double[] r, int par, int nThreads) {
        // vlastni vypocet, parallelni beh
        GmresThread[] parallel = new GmresThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new GmresThread(v, nThreads, x, r, par);
            parallel[v].start();
        }
        try {
            for (int v = 0; v < nThreads; v++) {
                parallel[v].join();
            }
        } catch (java.lang.InterruptedException e) {
            System.out.println(e);
        }
    }

    class GmresThread extends Thread {

        int nStart, nThreads, par, nt;
        double[] x, r;
        double residuumJacobi;

        GmresThread(int nStart, int nThreads, double[] x, double[] r, int par) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.x = x;
            this.r = r;
            this.par = par;
            nt = elems.length;
        }

        @Override
        public void run() {
            for (int i = nStart; i < nt; i = i + nThreads) {
                if (elems[i].insideComputeDomain) {
                    elems[i].residuumGmres(x, r, par);
                }
            }
        }
    }

    double ComputeResiduumJacobi(double[] x, double[] xn, int nThreads) {
        double residuum = 0;
        // vlastni vypocet, parallelni beh
        JacobiThread[] parallel = new JacobiThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new JacobiThread(v, nThreads, x, xn);
            parallel[v].start();
        }
        try {
            for (int v = 0; v < nThreads; v++) {
                parallel[v].join();
                residuum += parallel[v].residuumJacobi;
            }
        } catch (java.lang.InterruptedException e) {
            System.out.println(e);
        }

        return residuum;
    }

    class JacobiThread extends Thread {

        int nStart, nThreads, par, nt;
        double[] x, xn, r;
        double residuumJacobi;

        JacobiThread(int nStart, int nThreads, double[] x, double[] xn) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.x = x;
            this.xn = xn;
            nt = elems.length;
        }

        @Override
        public void run() {
            residuumJacobi = 0;
            for (int i = nStart; i < nt; i = i + nThreads) {
                if (elems[i].insideComputeDomain) {
                    residuumJacobi += elems[i].residuumJacobi(x, xn);
                }
            }
        }
    }

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
                
                case "MTJ":
                    //solver = new MTJsolver(elems, dofs, 500, par.iterativeSolverTol, par.nThreads);
                    break;    
                    
                default:
                    solver = new Gmres(elems, dofs, 50, 5, par.iterativeSolverTol, par.nThreads);
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
