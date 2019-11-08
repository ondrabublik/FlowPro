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

    abstract public void testDynamic(double dt, int newtonIter) throws IOException;

    public static MasterSolver factory(String simulationPath, Mesh[] meshes, Dynamics dyn,
            Equation eqn, Parameters par, State state, Domain domain, Object lock) {

        try {
            if (par.parallelMode) {
                switch (par.parallelSolverType.toLowerCase()) {
                    case "schwartz":
                        return new SchwartzImplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);

                    case "ksp":
                        return new KSPSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);
                }
            } else {
                switch (((meshes[0].getElems())[0].ti).getLocalSolverType()) {
                    case "localimplicit":
                        return new LocalImplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);

                    case "localexplicit":
                        return new LocalExplicitSolver(simulationPath, meshes, dyn, eqn, par, state, domain, lock);
                }
            }
        } catch (Exception e) {
            System.out.println("Unknown solver: " + ((meshes[0].getElems())[0].ti).getLocalSolverType());
        }

        return null;
    }
}
