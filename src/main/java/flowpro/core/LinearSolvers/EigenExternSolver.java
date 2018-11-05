package flowpro.core.LinearSolvers;

import eigenwrapper.EigenWrapper;
import flowpro.core.Parameters;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author obublik
 */
public class EigenExternSolver extends LinearSolver {

    EigenWrapper eigenWrap;

    EigenExternSolver(SparseMatrix A, Parameters par) {
        eigenWrap = new EigenWrapper();
    }

    @Override
    public boolean solve(double[] x, double[] b) {
        try {      
            eigenWrap.solveEigen(A.getRowIndexesCOO(), A.getColumnIndexesCOO(), A.getData(), b, x);
            
            return true;
        } catch (Exception ex) {
            Logger.getLogger(UmfpackExternSolver.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }
    }
}
