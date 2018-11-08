/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;

/**
 *
 * @author obublik
 */
public class ProstaIterace extends LinearSolver {
    
    int n, m, iterationLimit, nThreads;
    double tol;
    double[] r, aux;
    
    ProstaIterace(SparseMatrix A, Preconditioner M, int iterationLimit, double tol, int nThreads) {
        this.A = A;
        A.buildCRSformat();
        this.M = M;
        this.n = A.getDofs();
        this.iterationLimit = iterationLimit;
        this.tol = tol;
        this.nThreads = nThreads;
        
        r = new double[n];
        aux = new double[n];
    }
    
    public boolean solve(double[] x, double[] b){
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
    
    double norm(double[] a) {
        double n = 0;
        for (int i = 0; i < a.length; i++) {
            n = n + a[i] * a[i];
        }
        n = Math.sqrt(n);

        return n;
    }
}
