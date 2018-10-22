package flowpro.core.LinearSolvers;

import flowpro.core.Parameters;
import java.util.logging.Level;
import java.util.logging.Logger;
import umfpackwrapper.UmfpackWrapper;

/**
 *
 * @author obublik
 */
public class UmfpackExternSolver extends LinearSolver2 {

    UmfpackWrapper umfWrap;

    UmfpackExternSolver(SparseMatrix A, Parameters par) {
        A.buildCCSformat();
        umfWrap = new UmfpackWrapper();
    }

    @Override
    public boolean solve(double[] x, double[] b) {
        try {      
            umfWrap.solveUmfpack(A.getRowIndexesCCS(), A.getColumnIndexesCCS(), A.getData(), b, x);
            
            return true;
        } catch (Exception ex) {
            Logger.getLogger(ExternSolver.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }
    }
}
