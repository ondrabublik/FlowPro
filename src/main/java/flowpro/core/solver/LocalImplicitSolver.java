package flowpro.core.solver;

import flowpro.core.*;
import flowpro.api.Mat;
import flowpro.core.parallel.*;
import flowpro.api.Equation;
import static flowpro.core.FlowProMain.*;
import flowpro.core.element.Element;
import flowpro.api.Dynamics;
import flowpro.api.FluidForces;
import flowpro.core.LinearSolvers.LinearSolver;
import flowpro.core.element.Implicit;
import flowpro.core.meshDeformation.*;

import java.io.*;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.lang3.time.StopWatch;

/**
 * tento program resi Navierovy-Stokesovy rovnice na nestrukturovane siti pomoci
 * nespojite Galerkinovy metody
 */
public class LocalImplicitSolver extends MasterSolver {

    private static final Logger LOG = LoggerFactory.getLogger(MasterSolver.class);

//    private final Object lock;
    // common parameters
    private final Deformation dfm;
    private final Dynamics dyn;
    private final Equation eqn;
    private final Parameters par;

    // master and loner parameters
    private final State state;
    private final Object lock;
    private final String simulationPath;

    // slave and loner parameters
    private final Mesh mesh;
    private final Element[] elems;
    private final int dofs;
//    private Solution solution;

    /**
     * Constructor for master and local.
     *
     * @param simulationPath
     * @param meshes
     * @param dyn
     * @param eqn
     * @param par
     * @param state
     * @param domain
     * @param lock
     */
    public LocalImplicitSolver(String simulationPath, Mesh[] meshes, Dynamics dyn,
            Equation eqn, Parameters par, State state, Domain domain, Object lock) {
        this.simulationPath = simulationPath;
        this.dyn = dyn;
        this.eqn = eqn;
        this.par = par;
        this.state = state;
        mesh = meshes[0];
        elems = mesh.getElems();
        this.dofs = mesh.dofs;
        this.dfm = mesh.getDfm();
        this.lock = lock;
    }

    public Mesh getMesh() {
        return mesh;
    }

    private double timeStep(double CFL) {
        double dt = Double.MAX_VALUE;
        for (Element elem : elems) {
            double loc_dt = elem.delta_t(CFL);
            if (loc_dt < dt) {
                dt = loc_dt;
            }
        }

        return dt;
    }

    private double improveCFL(double cfl, double residuumRation) {
        if (par.varyCFL) {
            double beta = 1.2;
            return Math.max(Math.min(cfl * residuumRation, beta * cfl), cfl / beta);
        } else {
            return cfl;
        }
    }

    private void copyWo2W() {
        for (Element elem : elems) {
            elem.copyWo2W();
        }
    }

    private void copyW2Wo() {
        for (Element elem : elems) {
            elem.copyW2Wo();
        }
    }

    private double calculateResiduumW(double dt) {
        double resid = 0;
        for (Element elem : elems) {
            resid += elem.calculateResiduumW(dt);
        }

        return resid / elems.length;
    }

    private String infoToString(int totalSteps, double dt, long assembleTime, long solveTime) throws IOException {
        String timeStr = millisecsToTime(state.getOverallExecutionTime());

        return String.format("%d/%d  resid: %.2e,  dt: %.1e,  t: %.2f,  CFL: %1.2f,  CPU: %s, AT: %dms, ST: %dms",
                state.steps, totalSteps, state.residuum, dt, state.t, state.cfl, timeStr, assembleTime, solveTime);
    }

    public Solution solve() throws IOException {
        int nElems = elems.length;
        LinearSolver linSolver = LinearSolver.factory(elems, par);
        JacobiAssembler assembler = new JacobiAssembler(elems, par);
        double[] x = new double[dofs];
        StopWatch watch = new StopWatch();

        CFLSetup cflObj = new CFLSetup(par.cfl, par.varyCFL);
        double dto = -1;
        boolean converges = true;
        int totalSteps = state.steps + par.steps;
        long assembleTime = 0;
        long solveTime = 0;
        boolean isFirstIter = true;

        // save zero iteration
        if (par.animation && state.steps == 0) {
            Solution solution = mesh.getSolution();
            saveAnimationData(solution, state.steps);
        }

        LOG.info("computation has started...");
        watch.start();
        outerloop:
        for (++state.steps; state.steps <= totalSteps && state.residuum > par.residuum
                && state.t < par.tEnd; ++state.steps) {

            if (converges) {
                state.cfl = cflObj.getCFL(state.cfl, state.residuum);
            }
            // zastavovaci podminka
            if (state.cfl < (par.cfl / 20)) {
                LOG.error("algorithm does not converge - aborting computation");
                return mesh.getSolution();
            }

            // nastaveni dt
            double dt = timeStep(state.cfl);
            if (dto == -1) {
                dto = dt;
            }
            if (state.t + dt > par.tEnd) {
                dt = par.tEnd - state.t;
            }
            // set equation object state
            eqn.setState(state.t + dt, dt);

            //compute domain 
            if (par.solutionMonitorOn) {
                mesh.computeSolutionMonitor();
            }

            // compute artificial viscosity
            for (int i = 0; i < nElems; i++) {
                elems[i].limiter(isFirstIter);
            }

            assembleTime = 0;
            solveTime = 0;
            for (int s = 0; s < par.newtonIters; s++) {  // vnitrni iterace (Newton)
                // mesh deformation
                if (par.movingMesh) {
                    dfm.calculateForces(elems, dyn.getMeshMove());
                    dyn.computeBodyMove(dt, state.t, s, dfm.getFluidForces());
                    dfm.newMeshPositionAndVelocity(elems, elems[0].ti.getOrder(), dt, dto, dyn.getMeshMove());
                    if (isFirstIter) {
                        dfm.relaxFirstIteration(elems,dt);
                    }
                    dfm.recalculateMesh(elems, par.order);
                }

                long startTime = System.currentTimeMillis();
                assembler.assemble(dt, dto);
                assembleTime += System.currentTimeMillis() - startTime;

                // reseni soustavy rovnic
                Arrays.fill(x, 0.0);
                startTime = System.currentTimeMillis();
                converges = linSolver.solve(x);
                solveTime = +System.currentTimeMillis() - startTime;

                if (!converges) {
                    copyWo2W();
                    state.cfl = cflObj.reduceCFL(state.cfl);
                    --state.steps;
                    LOG.warn("GMRES does not converge, CFL reduced to " + state.cfl);
                    continue outerloop;
                }

                state.residuum = Mat.L1Norm(x) / nElems;
                if (par.newtonIters > 1) {
                    LOG.info("   " + s + "-inner iteration, reziduum = " + state.residuum);
                }

                // ulozeni novych hodnot
                for (int i = 0; i < nElems; i++) {
                    ((Implicit)elems[i].ti).updateW(x);
                }

                double iner_tol = 1e-4;  // zadat jako parametr !!!!!!!
                if (state.residuum < iner_tol) {
                    break;
                }
            }

            state.residuum = calculateResiduumW(dt);
            if (state.residuum == 0) {
                LOG.error(" computation error ");
                //break;
            }
            state.executionTime = watch.getTime();
            state.t += dt;
            dto = dt;
            saveResiduum(state.residuum, state.t, state.executionTime);
            
            copyW2Wo();
            mesh.updateTime(dt);
            if (par.movingMesh) {
                dfm.nextTimeLevel(elems);
                dyn.nextTimeLevel();
                dyn.savePositionsAndForces();
                for (int i = 0; i < nElems; i++) {
                    elems[i].nextTimeLevelMassMatrixes();
                }
            }

            String info = infoToString(totalSteps, dt, assembleTime, solveTime);
            if ((state.steps % par.saveRate) == 0) {
                Solution solution = mesh.getSolution();
                if (par.animation) {
                    saveAnimationData(solution, state.steps);
                }
                saveData(solution);
                LOG.info(info);
            }
            System.out.printf(info);
            System.out.println();
            isFirstIter = false;
        }
        watch.stop();

        --state.steps;
        state.executionTime = watch.getTime();
        LOG.info("computation has finished in " + millisecsToTime(state.executionTime));
        Solution solution = mesh.getSolution();

        return solution;
    }

    @Override
    public void testDynamic(double dt, int newtonIter) throws IOException {
        double t = 0;
        for (int step = 0; step <= par.steps; step++) {
            dyn.computeBodyMove(dt, t, newtonIter, new FluidForces(new double[2][dfm.nBodies], new double[1][dfm.nBodies], null, null, null));
            dyn.nextTimeLevel();
            dyn.savePositionsAndForces();
            t += dt;
            System.out.println(step + "-th iteration t = " + t);
        }
    }

    public void saveData(Solution sol) throws IOException {
        synchronized (lock) {
            state.save();
            eqn.saveReferenceValues(simulationPath + REF_VALUE_FILE_NAME);
            Mat.save(sol.avgW, simulationPath + "W.txt");
            Mat.save(sol.W, simulationPath + "We.txt");
            Mat.save(mesh.getArtificialViscosity(), simulationPath + "artificialViscosity.txt");
            if (par.movingMesh) {
                Mat.save(sol.vertices, simulationPath + "PXY.txt");
                Mat.save(sol.meshVelocity, simulationPath + "UXY.txt");
            }
            File[] content = new File(simulationPath + "output").listFiles();
            if (content != null) {
                for (File file : content) {
                    file.delete();
                }
            }

            lock.notify();
        }
        LOG.info("results have been saved into " + simulationPath);
    }

    public void saveAnimationData(Solution sol, int step) throws IOException {
        File directory = new File(simulationPath + "animation");
        if (!directory.exists()) {
            directory.mkdir();
        }
        synchronized (lock) {
            Mat.save(sol.avgW, simulationPath + "animation/W" + (10000000 + step) + ".txt");
            if (par.order > 1) {
                Mat.save(sol.W, simulationPath + "animation/We" + (10000000 + step) + ".txt");
            }
            if (par.movingMesh) {
                Mat.save(sol.vertices, simulationPath + "animation/vertices" + (10000000 + step) + ".txt");
                Mat.save(sol.meshVelocity, simulationPath + "animation/verticesVelocity" + (10000000 + step) + ".txt");
            }
            lock.notify();
        }
        LOG.info("results have been saved into " + simulationPath);
    }

    public class CFLSetup {

        double maxCFL;
        boolean varyCFL;
        double res = -1;
        double dres = 0;
        double alfa = 0.1;
        int logResOld;

        CFLSetup(double maxCFL, boolean varyCFL) {
            this.maxCFL = maxCFL;
            this.varyCFL = varyCFL;
        }

        double getCFL(double actualCFL, double residuum) {
            if (varyCFL) {
                if (res == -1) {
                    res = residuum;
                    logResOld = 1000;
                } else {
                    res = alfa * res + (1 - alfa) * residuum; // low pass filter
                }
                int logRes = (int) Math.log(res);
                if (logRes < logResOld) {
                    maxCFL *= 1.3;
                } else {
                    maxCFL *= 0.8;
                }
                logResOld = logRes;
                actualCFL += maxCFL / 5;
                if (actualCFL > maxCFL) {
                    actualCFL = maxCFL;
                }
                return actualCFL;
            } else {
                actualCFL += maxCFL / 20;
                if (actualCFL > maxCFL) {
                    actualCFL = maxCFL;
                }
                return actualCFL;
            }
        }

        double reduceCFL(double actualCFL) {
            return actualCFL / 1.5;
        }
    }

    public class JacobiAssembler {

        public Element[] elems;
        private final int nEqs;
        private final Parameters par;
        double dt, dt0;


        public JacobiAssembler(Element[] elems, Parameters par) {
            this.elems = elems;
            this.par = par;
            nEqs = elems[0].getNEqs();
        }

        // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
        public void assemble(double dt, double dto) {  // , int newtonIter
            this.dt = dt;
            this.dt0 = dt0;

            AssemblerThread[] assemblers = new AssemblerThread[par.nThreads];

            // vlastni vypocet, parallelni beh
            for (int v = 0; v < assemblers.length; v++) {
                assemblers[v] = new AssemblerThread(v);
                assemblers[v].start();
            }

            try {
                for (AssemblerThread assembler : assemblers) {
                    assembler.join();
                }
            } catch (java.lang.InterruptedException e) {
                System.err.println(e);
                System.exit(1);
            }
        }

        private class AssemblerThread extends Thread {

            private final int id;

            AssemblerThread(int id) {
                this.id = id;
            }

            @Override
            public void run() {
                for (int i = id; i < elems.length; i += par.nThreads) {
                    if (elems[i].insideComputeDomain) {
                        ((Implicit)elems[i].ti).assembleJacobiMatrix(dt,dt0);
                    }
                }
            }
        }
    }

    public void saveResiduum(double residuum, double t, double CPU) {
        try {
            FileWriter fw;
            fw = new FileWriter(simulationPath + "residuum.txt", true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter out = new PrintWriter(bw);
            out.println(Double.toString(residuum) + " " + Double.toString(t) + " " + Double.toString(CPU));
            out.close();
        } catch (Exception e) {
            //exception handling left as an exercise for the reader
        }
    }
}
