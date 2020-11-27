package flowpro.core.solver;

import flowpro.core.*;
import flowpro.api.Mat;
import flowpro.core.parallel.*;
import flowpro.api.Equation;
import static flowpro.core.FlowProMain.*;
import flowpro.api.Dynamics;
import flowpro.api.FluidForces;
import flowpro.core.DistributedLinearSolver.ParallelGmresMaster;
import flowpro.core.meshDeformation.*;
import flowpro.core.parallel.Domain.Subdomain;
import litempi.*;

import java.io.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.lang3.time.StopWatch;

/**
 * tento program resi Navierovy-Stokesovy rovnice na nestrukturovane siti pomoci
 * nespojite Galerkinovy metody
 */
public class KSPSolver extends MasterSolver {

    private static final Logger LOG = LoggerFactory.getLogger(MasterSolver.class);

//    private final Object lock;
    // common parameters
    private final Deformation dfm;
    private final Dynamics dyn;
    private final Equation eqn;
    private final Parameters par;

    // master and loner parameters
    private Mesh[] meshes;
    private final State state;
    private final Domain domain;
    private final Object lock;
    private final String simulationPath;
    private final MPIMaster mpi;

    // slave and loner parameters
    private final Mesh mesh;
//    private Solution solution;

    /**
     * Constructor for master and local.
     *
     * @param mpi
     * @param simulationPath
     * @param meshes
     * @param dyn
     * @param eqn
     * @param par
     * @param state
     * @param domain
     * @param lock
     */
    public KSPSolver(MPIMaster mpi, String simulationPath, Mesh[] meshes, Dynamics dyn,
            Equation eqn, Parameters par, State state, Domain domain, Object lock) {
        this.mpi = mpi;
        this.simulationPath = simulationPath;
        this.meshes = meshes;
        this.dyn = dyn;
        this.eqn = eqn;
        this.par = par;
        this.state = state;
        this.domain = domain;
        mesh = meshes[0];
        this.dfm = mesh.getDfm();
        this.lock = lock;
    }


    @Override
    public Mesh getMesh() {
        return mesh;
    }

    private String infoToString(int totalSteps, double dt, long assembleTime, long solveTime) throws IOException {
        String timeStr = millisecsToTime(state.getOverallExecutionTime());

        return String.format("%d/%d  resid: %.2e,  dt: %.1e,  t: %.2f,  CFL: %1.2f,  CPU: %s, AT: %dms, ST: %dms",
                state.steps, totalSteps, state.residuum, dt, state.t, state.cfl, timeStr, assembleTime, solveTime);
    }

    private Solution collectSolution(MPIMaster mpi) throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.SOLUTION_REQ));
        Solution[] sols = new Solution[mpi.nSlaves];
        for (int d = 0; d < mpi.nSlaves; ++d) {
            sols[d] = (Solution) mpi.receive(d, Tag.SOLUTION).getData();
        }
        return new Solution(sols, domain);
    }

    private void exchangeData(LiteElement[] liteElems, MPIMaster mpi) throws MPIException {
        // downloading data from the slave nodes into the central structure
        mpi.sendAll(new MPIMessage(Tag.DATA_REQ));
        for (int d = 0; d < domain.nDoms; ++d) {
            int[] mapL2G = domain.getSubdomain(d).mapL2G;
            LiteElement[] dataRcv = (LiteElement[]) mpi.receive(d, Tag.DATA_SLAVE_TO_MASTER).getData();
            for (LiteElement dataRcv1 : dataRcv) {
                //liteElems[mapL2G[dataRcv1.index]] = dataRcv1;
                liteElems[mapL2G[dataRcv1.index]] = new LiteElement(dataRcv1.index, dataRcv1.y);
            }
        }

        // uploading data from the central structure onto the slave nodes
        for (int d = 0; d < domain.nDoms; ++d) {
            Subdomain subdom = domain.getSubdomain(d);
            int[] mapG2L = subdom.mapG2L;
            int[] load = subdom.load;
            LiteElement[] dataSend = new LiteElement[load.length];
            for (int i = 0; i < load.length; i++) {
                dataSend[i] = new LiteElement(mapG2L[load[i]], liteElems[load[i]].y);
                //dataSend[i] = liteElems[load[i]];
                //dataSend[i].index = mapG2L[load[i]];
            }
            mpi.send(new MPIMessage(Tag.DATA_MASTER_TO_SLAVE, dataSend), d);
        }
    }

    @Override
    public Solution solve() throws MPIException, IOException {
        int nDoms = domain.nDoms;
 
        StopWatch watch = new StopWatch();
        double dto = 1;
        CFLSetup cflObj = new CFLSetup(par.cfl, par.varyCFL);
        long assembleTime = 0;
        long solveTime = 0;
        int tsr = 1; // for time save rate

        ParallelGmresMaster linSolver = new ParallelGmresMaster(par, 30, mpi, domain);

        try {
            LOG.info("sending initial data");
            for (int d = 0; d < nDoms; ++d) {
                mpi.send(new MPIMessage(Tag.INIT_DATA, meshes[d]), d);
                meshes[d] = null;  // garbage collector should clear the memory
            }
            meshes = null;  // garbage collector should clear the memory
            mpi.waitForAll(Tag.DATA_INITIALIZED);

            // central structure
            LiteElement[] liteElems = new LiteElement[domain.nElems];
            boolean convergesNewton = true;
            int totalSteps = state.steps + par.steps;

            LOG.info("computation has started...");

            // save zero iteration
            if (par.animation && state.steps == 0) {
                Solution solution = collectSolution(mpi);
                saveAnimationData(solution, state.steps);
            }

            watch.start();
            outerloop:
            for (++state.steps; state.steps <= totalSteps && state.residuum > par.residuum
                    && state.t < par.tEnd; ++state.steps) {
                if (convergesNewton) {
                    state.cfl = cflObj.getCFL(state.cfl, state.residuum);
                }
                // zastavovaci podminka
                if (state.cfl < (par.cfl / 20)) {
                    LOG.error("algorithm does not converge - aborting computation");
                    System.exit(1);
                }

                mpi.sendAll(new MPIMessage(Tag.TIME_STEP_REQ, state.cfl));
                double dt = Mat.min(mpi.receiveAllDouble(Tag.TIME_STEP));

                // solution monitor
                if (par.solutionMonitorOn) {
                    double[] monitorGlobal = null;
                    mpi.sendAll(new MPIMessage(Tag.GETSOLUTIONMONITOR));
                    for (int d = 0; d < nDoms; ++d) {
                        if (d == 0) {
                            monitorGlobal = (double[]) mpi.receive(d, Tag.LOCALSOLUTIONMONITOR).getData();
                        } else {
                            double[] monitorLocal = (double[]) mpi.receive(d, Tag.LOCALSOLUTIONMONITOR).getData();
                            for (int i = 0; i < monitorGlobal.length; i++) {
                                monitorGlobal[i] += monitorLocal[i];
                            }
                        }
                    }
                    mpi.sendAll(new MPIMessage(Tag.SETSOLUTIONMONITOR, (Object) monitorGlobal));
                    mpi.waitForAll(Tag.SOLUTIONMONITORSET);
                }

                for (int newtonIter = 0; newtonIter < par.newtonIters; newtonIter++) {
                    // mesh deformation
                    if (par.movingMesh) {
                        mpi.sendAll(new MPIMessage(Tag.ALE_CALCULATE_FORCES, dyn.getMeshMove()));
                        ForcesAndDisplacements[] forDis = new ForcesAndDisplacements[nDoms];
                        for (int d = 0; d < nDoms; ++d) {
                            forDis[d] = (ForcesAndDisplacements) mpi.receive(d, Tag.FORCES).getData();
                        }
                        dyn.computeBodyMove(dt, state.t, newtonIter, forDis[0].combine(forDis));
                        mpi.sendAll(new MPIMessage(Tag.ALE_NEW_MESH_POSITION, new ForcesAndDisplacements(dt, dto, dyn.getMeshMove())));
                        mpi.waitForAll(Tag.MESH_POSITION_UPDATED);
                    }

                    // assembling
                    long startTime = System.currentTimeMillis();
                    mpi.sendAll(new MPIMessage(Tag.ASSEMBLE, dt));
                    mpi.waitForAll(Tag.ASSEMBELED);
                    assembleTime = System.currentTimeMillis() - startTime;

                    // solve
                    startTime = System.currentTimeMillis();
                    convergesNewton = linSolver.solve();
                    solveTime = System.currentTimeMillis() - startTime;

                    // data transfer
                    exchangeData(liteElems, mpi);
                    mpi.waitForAll(Tag.DATA_UPDATED);

                    if (!convergesNewton) {
                        mpi.sendAll(new MPIMessage(Tag.PREVIOUS_TIME_LEVEL));
                        state.cfl = cflObj.reduceCFL(state.cfl);
                        --state.steps;
                        LOG.warn("GMRES does not converge, CFL reduced to " + state.cfl);
                        continue outerloop;
                    }

                    mpi.sendAll(new MPIMessage(Tag.UPDATE_NEWTON));
                    mpi.waitForAll(Tag.NEWTON_UPDATED);
                }

                mpi.sendAll(new MPIMessage(Tag.LIMITER_and_NEXT_TIME_LEVEL));
                double resid = 0.0;
                for (int d = 0; d < nDoms; ++d) {
                    resid += (double) mpi.receive(d, Tag.RESIDUUM).getData();
                }
                state.residuum = resid/domain.nElems;
                if (state.residuum == 0) {
                    LOG.error("computation error");
                    break;
                }

                if (par.movingMesh) {
                    mpi.sendAll(new MPIMessage(Tag.ALE_NEXT_TIME_LEVEL));
                    mpi.waitForAll(Tag.NEW_MESH_POSITION_UPDATED);
                    dyn.nextTimeLevel();
                    dyn.savePositionsAndForces();
                }

                dto = dt;
                state.t += dt;
                state.executionTime = watch.getTime();
                saveResiduum(state.residuum, state.t, state.executionTime);
                String info = infoToString(totalSteps, dt, assembleTime, solveTime);
                if ((state.steps % par.saveRate) == 0 || (state.t > tsr*par.timeSaveRate)) {
                    Solution solution = collectSolution(mpi);
                    if (par.animation) {
                        saveAnimationData(solution, state.steps);
                    }
                    saveData(solution);
                    LOG.info(info);
                    tsr++;
                }
                System.out.printf(info);
                System.out.println();
                if ((state.steps % 50) == 0) {
                    mpi.sendAll(new MPIMessage(Tag.RESET_OUTPUT_WRITER));
                    mpi.reset();
                }
            }
            watch.stop();
            state.executionTime = watch.getTime();
            --state.steps;
            LOG.info("computation has finished in " + millisecsToTime(state.executionTime));
            LOG.info("collecting data");
            Solution solution = collectSolution(mpi);

            return solution;
        } finally {
            mpi.sendAll(new MPIMessage(Tag.CLOSE));
            mpi.close();
        }
    }

    @Override
    public void testDynamic(double dt, int newtonIter) throws IOException {
        double t = 0;
        for (int step = 0; step <= par.steps; step++) {
			FluidForces[] fluidForces = new FluidForces[dfm.nBodies];
			for (int b = 0; b < dfm.nBodies; b++) {
				fluidForces[b] = new FluidForces(new double[2], new double[1]);
			}
            dyn.computeBodyMove(dt, t, newtonIter, fluidForces);			
            dyn.nextTimeLevel();
            dyn.savePositionsAndForces();
            t += dt;
            System.out.println(step + "-th iteration t = " + t);
        }
    }

    @Override
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
            Mat.save(sol.avgW, simulationPath + "animation/W" + String.format("%08d", step) + ".txt");
            if (par.order > 1) {
                Mat.save(sol.W, simulationPath + "animation/We" + String.format("%08d", step) + ".txt");
            }
            if (par.movingMesh) {
                Mat.save(sol.vertices, simulationPath + "animation/vertices" + String.format("%08d", step) + ".txt");
                Mat.save(sol.meshVelocity, simulationPath + "animation/verticesVelocity" + String.format("%08d", step) + ".txt");
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

        public CFLSetup(double maxCFL, boolean varyCFL) {
            this.maxCFL = maxCFL;
            this.varyCFL = varyCFL;
        }

        public double getCFL(double actualCFL, double residuum) {
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

    public void saveResiduum(double residuum, double t, double CPU) {
        try {
            FileWriter fw;
            fw = new FileWriter(simulationPath + "residuum.txt", true);
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter out = new PrintWriter(bw);
            out.println(Double.toString(residuum) + " " + Double.toString(t) + " " + Double.toString(CPU));
            out.close();
        } catch (IOException ex) {
            LOG.warn("cannot write into file {}: {}", "residuum.txt", ex.getMessage());
        }
    }
}

