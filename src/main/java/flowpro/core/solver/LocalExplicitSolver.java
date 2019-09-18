package flowpro.core.solver;

import flowpro.core.*;
import flowpro.api.Mat;
import flowpro.core.parallel.*;
import flowpro.api.Equation;
import static flowpro.core.FlowProMain.*;
import flowpro.core.element.Element;
import flowpro.api.Dynamics;
import flowpro.api.FluidForces;
import flowpro.core.element.RK3Element;
import flowpro.core.meshDeformation.*;

import java.io.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.lang3.time.StopWatch;

/**
 * tento program resi Navierovy-Stokesovy rovnice na nestrukturovane siti pomoci
 * nespojite Galerkinovy metody
 */
public class LocalExplicitSolver extends MasterSolver {

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
    public LocalExplicitSolver(String simulationPath, Mesh[] meshes, Dynamics dyn,
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

    private void copyW2Wo() {
        for (Element elem : elems) {
            elem.copyW2Wo();
        }
    }

    private void computeInvertMassMatrix() {
        for (Element elem : elems) {
            elem.computeInvertMassMatrix();
        }
    }

    private double calculateResiduumW(double dt) {
        double resid = 0;
        for (Element elem : elems) {
            resid += elem.calculateResiduumW(dt);
        }

        return resid / elems.length;
    }

    private String infoToString(int totalSteps, double dt) throws IOException {
        String timeStr = millisecsToTime(state.getOverallExecutionTime());

        return String.format("%d/%d  resid: %.2e,  dt: %.1e,  t: %.2f,  CFL: %1.2f,  CPU: %s",
                state.steps, totalSteps, state.residuum, dt, state.t, state.cfl, timeStr);
    }

    public Solution solve() throws IOException {
        if (par.movingMesh) {
            LOG.error("Moving mesh not supported for local time stepping explicit method!");
        }

        if (par.orderInTime == 1) {
            LOG.info("Explicit first order time integration method is not suported! Setting time order to second order!");
            par.orderInTime = 2;
        }

        computeInvertMassMatrix();
        LocalTimeStepIterator ltsIter = new LocalTimeStepIterator(elems, par, state.t);
        StopWatch watch = new StopWatch();

        CFLSetup cflObj = new CFLSetup(par.cfl, par.varyCFL);
        int totalSteps = state.steps + par.steps;

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

            state.cfl = cflObj.getCFL(state.cfl, state.residuum);

            // nastaveni dt
            double dt = timeStep(state.cfl);
            if (state.t + dt > par.tEnd) {
                dt = par.tEnd - state.t;
            }

            //compute domain 
            if (par.solutionMonitorOn) {
                mesh.computeSolutionMonitor();
            }

            // computation
            ltsIter.iterate(1, dt);
            ltsIter.iterate(11, dt);
            ltsIter.iterate(2, dt);
            ltsIter.iterate(21, dt);
            ltsIter.iterate(3, dt);
            ltsIter.iterate(31, dt);

            state.residuum = calculateResiduumW(dt);
            state.executionTime = watch.getTime();
            state.t += dt;
            saveResiduum(state.residuum, state.t, state.executionTime);
            copyW2Wo();

            String info = infoToString(totalSteps, dt);
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
        }
        watch.stop();

        --state.steps;
        state.executionTime = watch.getTime();
        LOG.info("computation has finished in " + millisecsToTime(state.executionTime));
        Solution solution = mesh.getSolution();

        return solution;
    }

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
            }
            lock.notify();
        }
        LOG.info("results have been saved into " + simulationPath);
    }

    class CFLSetup {

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

    class LocalTimeStepIterator {

        public Element[] elems;
        private final int nEqs;
        private final Parameters par;

        LocalTimeStepIterator(Element[] elems, Parameters par, double t) {
            this.elems = elems;
            this.par = par;
            nEqs = elems[0].getNEqs();
        }

        // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
        public void iterate(int step, double dt) {  // , int newtonIter

            LocalTimeStepIteratorIteratorThread[] iteratorThrds = new LocalTimeStepIteratorIteratorThread[par.nThreads];

            // vlastni vypocet, parallelni beh
            for (int v = 0; v < iteratorThrds.length; v++) {
                iteratorThrds[v] = new LocalTimeStepIteratorIteratorThread(v, step, dt);
                iteratorThrds[v].start();
            }

            try {
                for (LocalTimeStepIteratorIteratorThread iter : iteratorThrds) {
                    iter.join();
                }
            } catch (java.lang.InterruptedException e) {
                System.err.println(e);
                System.exit(1);
            }
        }
    }

    private class LocalTimeStepIteratorIteratorThread extends Thread {

        private final int id;
        private final int step;
        private final double dt;

        LocalTimeStepIteratorIteratorThread(int id, int step, double dt) {
            this.id = id;
            this.step = step;
            this.dt = dt;
        }

        @Override
        public void run() {
            for (int i = id; i < elems.length; i += par.nThreads) {
                if (elems[i].insideComputeDomain) {
                    ((RK3Element) elems[i].ti).computeExplicitStep(step, dt);
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
