package flowpro.core.LinearSolvers;

import flowpro.core.Mesh.*;

/**
 *
 * @author obublik
 */
public class Jacobi extends LinearSolver {

    int n;
    int iterationLimit;
    int nThreads;
    double tol;

    double[] xn;

    public Jacobi(Element[] elems, int n, int iterationLimit, double tol, int nThreads) {
        this.elems = elems;
        this.n = n;
        this.iterationLimit = iterationLimit;
        this.tol = tol;
        this.nThreads = nThreads;

        xn = new double[n];
    }

    @Override
    public boolean solve(double[] x) {
        for (int iter = 0; iter < iterationLimit; iter++) {
            double error = Math.sqrt(ComputeResiduumJacobi(x, xn, nThreads));
            System.arraycopy(xn, 0, x, 0, n);
            
            if(Math.sqrt(error) < tol){
                return true;
            }
        }

        return false;
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
}