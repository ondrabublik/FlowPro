/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.Parameters;
import litempi.*;

/**
 *
 * @author obublik
 */
public class ParallelGmresMaster extends LinearSolver{

    int iterationLimit, nThreads;
    double tol;

    public ParallelGmresMaster(Parameters par, MPIMaster mpi) {
        iterationLimit = 10;
        tol = par.iterativeSolverTol;
        nThreads = par.nThreads;
    }

    public boolean solve(double[] x, double[] b) {
        boolean converged = false;
        for(int i = 0; i < iterationLimit; i++){
            multiplyAndUpdate();
            dataExchange();
            double error = norm();
            if(error < tol){
                converged = true;
                break;
            }
        }
        return converged;
    }
}
