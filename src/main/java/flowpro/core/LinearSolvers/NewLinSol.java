/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;
import flowpro.core.Mesh;
import flowpro.core.Parameters;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class NewLinSol extends LinearSolver {

    int dofs;
    double[] b;

    SparseMatrix A;
    Preconditioner M;
    Gmres2 solver;

    NewLinSol(Mesh.Element[] elems, int dofs, Parameters par) throws IOException {
        this.elems = elems;
        this.dofs = dofs;
        b = new double[dofs];
        A = new SparseMatrix(elems);
        A.buildCRSformat();
        M = Preconditioner.factory(par);
        M.setMatrix(A);
        solver = new Gmres2(A, M, 30, 5, par.iterativeSolverTol, par.nThreads);
    }

    @Override
    public boolean solve(double[] x) {
        
        A.updateData(); // update data in matrix A
        A.updateB(b);   // update RHS
        M.factor();     // update preconditioner M
        
        boolean flag = solver.solve(x, b);
        return flag;
    }
}
