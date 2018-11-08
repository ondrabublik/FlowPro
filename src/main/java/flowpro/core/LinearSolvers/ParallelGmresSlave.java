package flowpro.core.LinearSolvers;

import static flowpro.core.LinearSolvers.LinearSolver.A;
import static flowpro.core.LinearSolvers.LinearSolver.M;
import static flowpro.core.LinearSolvers.LinearSolver.b;
import flowpro.core.LinearSolvers.preconditioners.Preconditioner;
import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import litempi.MPIMessage;

/**
 *
 * @author obublik
 */
public class ParallelGmresSlave extends LinearSolver {

    int n, m, iterationLimit, nThreads;
    double tol;
    double[][] V, H;
    double[] cs, sn, e1, w, r, aux;

    //public ParallelGmresSlave(SparseMatrix A, Preconditioner M, int m, int iterationLimit, double tol, int nThreads) {
    public ParallelGmresSlave(Element[] elems, Parameters par) {
        // build matrix structure
        A = new SparseMatrix(elems);
        A.buildCRSformat();

        // define preconditiner
        M = Preconditioner.factory(par);
        M.setMatrix(A);

        // alocate RHS
        b = new double[A.getDofs()];

        n = A.getDofs();
        m = 30;
        iterationLimit = 10;
        tol = par.iterativeSolverTol;
        nThreads = par.nThreads;

        // initialize workspace
        V = new double[m + 1][n];
        H = new double[m + 1][m];
        cs = new double[m];
        sn = new double[m];
        e1 = new double[m + 1];
        w = new double[n];
        r = new double[n];
        aux = new double[n];
    }

    public MPIMessage doWork(MPIMessage mpi) {
        return null;
    }

    public boolean solve(double[] x, double[] b) {
        boolean converged = false;
        for(int i = 0; i < iterationLimit; i++){
            A.SubstrMult(aux, b, x, nThreads);
            M.apply(r, aux);
            for(int j = 0; j < n; j++){
                x[j] += r[j];
            }
            double error = norm(aux);
            if(error < tol){
                converged = true;
                break;
            }
        }
        return converged;
    }

    void updateSolution(double[] y, double[] x, int i, int nThreads) {
        // vlastni vypocet, parallelni beh
        UpdateSolutionThread[] parallel = new UpdateSolutionThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new UpdateSolutionThread(v, nThreads, y, x, i);
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

    class UpdateSolutionThread extends Thread {

        int nStart, nThreads, i, n;
        double[] y, x;

        UpdateSolutionThread(int nStart, int nThreads, double[] y, double[] x, int i) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.y = y;
            this.x = x;
            this.i = i;
            n = x.length;
        }

        @Override
        public void run() {
            // paralelni sestavovani a plneni matic
            for (int j = nStart; j < n; j = j + nThreads) {
                for (int k = 0; k <= i; k++) {
                    x[j] = x[j] + V[k][j] * y[k];
                }
            }
        }
    }

    double[] rotmat(double a, double b) {
        // Compute the Givens rotation matrix parameters for a and b.
        double c, s, temp;
        if (b == 0) {
            c = 1;
            s = 0;
        } else if (Math.abs(b) > Math.abs(a)) {
            temp = a / b;
            s = 1 / Math.sqrt(1 + temp * temp);
            c = temp * s;
        } else {
            temp = b / a;
            c = 1 / Math.sqrt(1 + temp * temp);
            s = temp * c;
        }
        return new double[]{c, s};
    }

    double[] vectorScalarProduct(double[] a, double b) {
        double[] c = new double[n];
        for (int i = 0; i < a.length; i++) {
            c[i] = a[i] * b;
        }
        return c;
    }

    double norm(double[] a) {
        double n = 0;
        for (int i = 0; i < a.length; i++) {
            n = n + a[i] * a[i];
        }
        n = Math.sqrt(n);

        return n;
    }

    double scalarProduct(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            s = s + a[i] * b[i];
        }
        return s;
    }

    // Gaussian elimination with partial pivoting
    double[] lsolve(double[][] A, double[] b, int N) {
        // int N  = b.length;
        for (int p = 0; p < N; p++) {
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }
}
