/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.solver;

import flowpro.api.Dynamics;
import flowpro.api.Equation;
import flowpro.core.Mesh;
import flowpro.core.Parameters;
import flowpro.core.Solution;
import flowpro.core.State;
import flowpro.core.parallel.Domain;
import java.io.IOException;
import litempi.MPIException;

/**
 *
 * @author obublik
 */
abstract public class MasterSolver {

    abstract public Solution solve() throws MPIException, IOException;
    abstract public void saveData(Solution sol) throws IOException;
    abstract public Mesh getMesh();
    abstract public void testDynamic(double dt) throws IOException;
    
    public static MasterSolver factory(String simulationPath, Mesh[] meshes, Dynamics dyn,
            Equation eqn, Parameters par, State state, Domain domain, Object lock) {

        try {
            switch (par.solverType.toLowerCase()) {
                case "localimplicit":
                    return new LocalImplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);

                case "localexplicit":
                    par.isExplicit = true;
                    return new LocalExplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);
                
                case "distimplicitschwartz":
                    return new SchwartzImplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);
                    
                case "distKSP":
                    return new KSPSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);
                
            }
        } catch (Exception e) {
            System.out.println("Unknown solver: " + par.solverType);
        }
        return null;
    }
}
