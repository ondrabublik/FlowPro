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
    
    public static SlaveSolver factory(String parallelSolverType, String masterIP, int masterPort) {

        try {
            switch (parallelSolverType.toLowerCase()) {              
                case "distschwartz":
                    return new SchwartzImplicitSolverSlave(masterIP, masterPort);
                case "distksp":
                    return new KSPSolverSlave(masterIP, masterPort);
                
            }
        } catch (Exception e) {
            System.out.println("Unknown slave solver: " + parallelSolverType);
        }
        return null;
    }
}
