/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.solver;

import java.io.IOException;
import litempi.MPIException;

/**
 *
 * @author obublik
 */
abstract public class SlaveSolver {

    abstract public void solve() throws MPIException, IOException;
    
    public static SlaveSolver factory(String parallelSolverType, int slavePort) throws IOException, MPIException {
        switch (parallelSolverType.toLowerCase()) {
            case "schwartz":
                return new SchwartzImplicitSolverSlave(slavePort);
            case "ksp":
                return new KSPSolverSlave(slavePort);
            default:
                throw new IOException("unknown slave solver " + parallelSolverType);
        }
    }
}
