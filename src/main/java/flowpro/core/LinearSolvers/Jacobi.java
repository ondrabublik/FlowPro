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
}