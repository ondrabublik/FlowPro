package flowpro.core;

import flowpro.api.Mat;
import flowpro.core.parallel.*;
import flowpro.core.LinearSolvers.LinearSolver;
import flowpro.api.Equation;
import static flowpro.core.FlowProMain.*;
import flowpro.core.Mesh.Element;
import fetcher.FetcherServer;
import fetcher.ZipFile;
import flowpro.api.Dynamics;
import flowpro.api.FluidForces;
import flowpro.api.MeshMove;
import flowpro.core.meshDeformation.*;
import flowpro.core.parallel.Domain.Subdomain;
import flowpro.core.timeIntegration.TimeIntegration;
import litempi.*;

import java.io.*;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.lang3.time.StopWatch;

/**
 * tento program resi Navierovy-Stokesovy rovnice na nestrukturovane siti pomoci
 * nespojite Galerkinovy metody
 */
public class Solver {

    private static final Logger LOG = LoggerFactory.getLogger(Solver.class);
    private static final StopWatch tempWatch = new StopWatch();
    private static final StopWatch tempWatch2 = new StopWatch();

//    private final Object lock;
    // common parameters
    private final Deformation dfm;
    private final Dynamics dyn;
    private final Equation eqn;
    private final TimeIntegration tim;
    private final Parameters par;

    // master and loner parameters
    private Mesh[] meshes;
    private State state;
    private Domain domain;
    private final Object lock;
    private String simulationPath;

    // slave and loner parameters
    private Mesh mesh;
    private Element[] elems;
    private MPISlave mpiSlave;
    private int dofs;
//    private Solution solution;

    /**
     * Constructor for master and loner.
     *
     * @param simulationPath
     * @param meshes
     * @param dyn
     * @param eqn
     * @param tim
     * @param par
     * @param state
     * @param domain
     * @param lock
     */
    public Solver(String simulationPath, Mesh[] meshes, Dynamics dyn,
            Equation eqn, TimeIntegration tim, Parameters par, State state, Domain domain, Object lock) {
        this.simulationPath = simulationPath;
        this.meshes = meshes;
        this.dyn = dyn;
        this.eqn = eqn;
        this.tim = tim;
        this.par = par;
        this.state = state;
        this.domain = domain;
        mesh = meshes[0];
        elems = mesh.getElems();
        this.dofs = mesh.dofs;
        this.dfm = mesh.getDfm();
        this.lock = lock;
    }

    /**
     * Constructor for slave.
     *
     * @param masterIP
     * @param masterPort
     * @throws IOException
     * @throws MPIException
     */
    public Solver(String masterIP, int masterPort) throws IOException, MPIException {
        mpiSlave = new MPISlave(masterIP, masterPort);
        MPIMessage msg = mpiSlave.receive();

        if (msg.tag != Tag.INIT_DATA) {
            throw new MPIException("expected InitMessage");
        }

        mesh = (Mesh) msg.getData();
        meshes = new Mesh[]{mesh};
        mesh.init();
        elems = mesh.getElems();

        par = mesh.getPar();
        eqn = mesh.getEqn();
        tim = mesh.getTim();
        dfm = mesh.getDfm();
        dyn = null;
        dofs = mesh.dofs;
        lock = null;

        mpiSlave.send(new MPIMessage(Tag.DATA_INITIALIZED));
    }

    public Mesh getMesh() {
        return mesh;
    }

    private double timeStep(double CFL) {
        double dt = Double.MAX_VALUE;
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                double loc_dt = elem.delta_t(CFL);
                if (loc_dt < dt) {
                    dt = loc_dt;
                }
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

    private void computeArtificialViscosity() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.limiter();
            }
        }
    }

    private void callLimiters() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.limitUnphysicalValues();
            }
        }
    }

    private void getTimeLevelBack() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.timeStorage.getTimeLevelBack(elem.W);
            }
        }
    }

    private void saveTimeLevel(double t) {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.timeStorage.saveTimeLevel(t, elem.M, elem.W, elem.RHS_loc);
            }
        }
    }

    private void computeInvertMassMatrix() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.computeInvertMassMatrix();
            }
        }
    }

    private void updateRHS(double x[]) {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.updateRHS(x);
            }
        }
    }

    private void updateW(double x[]) {
        for (Element elem : elems) {
            elem.updateW(x);
        }
    }

    private double calculateResiduumW(double dt) {
        double resid = 0;
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                resid += elem.timeStorage.calculateResiduumW(dt, elem.W);
            }
        }
        return resid / elems.length;
    }

    public void slaveSolve() throws IOException, MPIException {
        LOG.info("slave is running...");
        MPISlave mpi = mpiSlave;
        mesh = meshes[0];
        try {
            LinearSolver linSolver = LinearSolver.factory(par, elems, dofs);
            JacobiAssembler assembler = new JacobiAssembler(elems, par);
            double[] x = new double[dofs];
            double[] y = new double[dofs];
            double dt = -1.0;
            double dto;
            boolean firstTimeStep = true;

            while (true) {
                /* receive message */
                MPIMessage inMsg = mpi.receive();
//                LOG.debug("received {}", inMsg.tag);
                MPIMessage outMsg = null;

                switch (inMsg.tag) {
                    case Tag.TIME_STEP_REQ: // vypocet dt
                        double CFL = (double) inMsg.getData();
                        double locDt = timeStep(CFL);
                        computeArtificialViscosity();

                        outMsg = new MPIMessage(Tag.TIME_STEP, locDt);
                        break;

                    case Tag.GETSOLUTIONMONITOR:
                        //compute solution monitor
                        mesh.computeSolutionMonitor();
                        outMsg = new MPIMessage(Tag.LOCALSOLUTIONMONITOR, mesh.integralMonitor);
                        break;

                    case Tag.SETSOLUTIONMONITOR:
                        //compute solution monitor
                        mesh.integralMonitor = (double[]) inMsg.getData();
                        mesh.solMonitor.combineMonitoredValues(mesh.integralMonitor);
                        outMsg = new MPIMessage(Tag.SOLUTIONMONITORSET);
                        break;

                    case Tag.ALE_CALCULATE_FORCES:
                        MeshMove[] mshMov = (MeshMove[]) inMsg.getData();
                        dfm.calculateForces(elems, mshMov);
                        outMsg = new MPIMessage(Tag.FORCES, new ForcesAndDisplacements(dfm.getFluidForces()));
                        break;

                    case Tag.ALE_NEW_MESH_POSITION:
                        ForcesAndDisplacements disp = (ForcesAndDisplacements) inMsg.getData();
                        dfm.newMeshPosition(elems, par.orderInTime, disp.getDto(), disp.getDt(), disp.getMeshMove());
                        if (firstTimeStep) {
                            dfm.relaxFirstIteration(elems);
                        }
                        dfm.recalculateMesh(elems, par.order);
                        outMsg = new MPIMessage(Tag.MESH_POSITION_UPDATED);
                        break;

                    case Tag.ALE_NEXT_TIME_LEVEL:
                        dfm.nextTimeLevel(elems);
                        outMsg = new MPIMessage(Tag.NEW_MESH_POSITION_UPDATED);
                        break;

                    case Tag.ASSEMBLE_AND_SOLVE: // inicialization if equation Ax=b
                        // dalo by se tomu vyhnout kdybysme ukladali souboru dt a dto
                        if (firstTimeStep) {  //  if dt has not yet been set
                            dt = (double) inMsg.getData();
                            dto = dt;
                            firstTimeStep = false;
                        } else {  // if dt has already been set
                            dto = dt;
                            dt = (double) inMsg.getData();
                        }
                        assembler.assemble(dt, dto);
                        Arrays.fill(y, 0.0);
                    // NO BREAK! continue to the following tag

                    case Tag.SOLVE: // reseni rovnice Ax = b - Ay
                        Arrays.fill(x, 0.0);
                        boolean converges = linSolver.solve(x);
                        Mat.plusEqual(y, x);
                        updateRHS(x);  // b = b - Ax
                        double schwarzResid = Mat.L1Norm(x) / x.length;

                        outMsg = new MPIMessage(Tag.GMRES_RESULT, new GmresResult(converges, schwarzResid));
                        break;

                    case Tag.UPDATE_NEWTON: // update newton
                        updateW(y);
                        outMsg = new MPIMessage(Tag.NEWTON_UPDATED);
                        break;

                    case Tag.LIMITER_and_NEXT_TIME_LEVEL: // next time level
                        double t = (double) inMsg.getData();
                        double residuum = calculateResiduumW(dt);
                        callLimiters();
                        saveTimeLevel(t);
                        mesh.updateTime(dt);
                        outMsg = new MPIMessage(Tag.RESIDUUM, residuum);
                        break;

                    case Tag.DATA_REQ: // saving data into central structure
                        int nSave = mesh.save.length;
                        LiteElement[] dataSend = new LiteElement[nSave];
                        for (int i = 0; i < nSave; i++) {
                            int j = mesh.save[i];
                            double[] yElem = new double[mesh.nEqs * elems[j].nBasis];
                            for (int m = 0; m < mesh.nEqs; m++) {
                                for (int p = 0; p < elems[j].nBasis; p++) {
                                    int ind = elems[j].nBasis * m + p;
                                    yElem[ind] = y[elems[j].gi_U[ind]];
                                }
                            }
                            dataSend[i] = new LiteElement(j, yElem);
                        }
                        outMsg = new MPIMessage(Tag.DATA_SLAVE_TO_MASTER, dataSend);
                        break;

                    case Tag.DATA_MASTER_TO_SLAVE: // loading data from central structure
                        LiteElement[] dataReceive = (LiteElement[]) inMsg.getData();
                        for (LiteElement dataReceive1 : dataReceive) {
                            int j = dataReceive1.index;
                            for (int m = 0; m < mesh.nEqs; m++) {
                                for (int p = 0; p < elems[j].nBasis; p++) {
                                    int ind = elems[j].nBasis * m + p;
                                    y[elems[j].gi_U[ind]] = dataReceive1.y[ind];
                                }
                            }
                        }
                        outMsg = new MPIMessage(Tag.DATA_UPDATED);
                        break;

                    case Tag.SOLUTION_REQ: // zaverecne odeslani vsech dat (krome vnejsiho prekrivu)
                        Solution localSolution = mesh.getSolution();
//                        Solution localSolution = new Solution(mesh.getW(), mesh.getAvgW(), mesh.getDetailW());
                        LOG.info("sending solution");
                        outMsg = new MPIMessage(Tag.SOLUTION, localSolution);
                        break;

                    case Tag.PREVIOUS_TIME_LEVEL: // GMRES does not converge, go to previous time level
                        getTimeLevelBack();
                        break;

                    case Tag.RESET_OUTPUT_WRITER:
                        mpi.reset();
                        break;

                    case Tag.CLOSE: // Ukonceni
                        LOG.info("received command to shut down");
                        return;
                }

                /* send message */
                if (outMsg != null) {
//                    LOG.debug("sending {}", outMsg.tag);
                    mpi.send(outMsg);
                }
            }
        } catch (MPIException ex) {
            throw ex;
        } catch (Throwable th) {
            mpi.send(new MPIExceptionMessage(th));
            throw th;
        } finally {
            mpi.close();
        }
    }

    private String infoToString(int totalSteps, double dt) throws IOException {
        String timeStr = millisecsToTime(state.getOverallExecutionTime());

        return String.format("%d/%d  resid: %.2e,  dt: %.1e,  t: %.2f,  CFL: %1.2f,  CPU: %s",
                state.steps, totalSteps, state.residuum, dt, state.t, state.cfl, timeStr);
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
            tempWatch.resume();
            LiteElement[] dataRcv = (LiteElement[]) mpi.receive(d, Tag.DATA_SLAVE_TO_MASTER).getData();
            tempWatch.suspend();
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
            tempWatch2.resume();
            mpi.send(new MPIMessage(Tag.DATA_MASTER_TO_SLAVE, dataSend), d);
            tempWatch2.suspend();
        }
    }

    public Solution masterSolve() throws MPIException, IOException {
        int nDoms = domain.nDoms;
        runFetcher(nDoms, par.fetcherPort, par.masterIP, par.masterPort);
        MPIMaster mpi = new MPIMaster(nDoms, par.masterIP, par.masterPort);
        StopWatch watch = new StopWatch();
        StopWatch transferWatch = new StopWatch();
        tempWatch.start();
        tempWatch.suspend();
        tempWatch2.start();
        tempWatch2.suspend();
        double dto = 1;
        CFLSetup cflObj = new CFLSetup(par.cfl, par.varyCFL);

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
            boolean convergesSchwarz = true;
            int totalSteps = state.steps + par.steps;

            LOG.info("computation has started...");

            // save zero iteration
            if (par.animation && state.steps == 0) {
                Solution solution = collectSolution(mpi);
                saveAnimationData(solution, state.steps);
            }

            watch.start();
            transferWatch.start();
            transferWatch.suspend();
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
                        dyn.computeBodyMove(dt, state.t, forDis[0].combine(forDis));
                        mpi.sendAll(new MPIMessage(Tag.ALE_NEW_MESH_POSITION, new ForcesAndDisplacements(dt, dto, dyn.getMeshMove())));
                        mpi.waitForAll(Tag.MESH_POSITION_UPDATED);
                    }

                    mpi.sendAll(new MPIMessage(Tag.ASSEMBLE_AND_SOLVE, dt));
                    convergesNewton = true;
                    for (int d = 0; d < nDoms; ++d) {
                        GmresResult results = (GmresResult) mpi.receive(d, Tag.GMRES_RESULT).getData();
                        if (results.converges == false) {
                            convergesNewton = false;
                        }
                    }

                    transferWatch.resume();
                    exchangeData(liteElems, mpi);
                    mpi.waitForAll(Tag.DATA_UPDATED);
                    transferWatch.suspend();

                    if (convergesNewton) {
                        double schwarzResid = Double.MAX_VALUE;
                        for (int schwarzIter = 1; schwarzIter < par.schwarzIters
                                && schwarzResid > par.schwarzTol; ++schwarzIter) {

                            mpi.sendAll(new MPIMessage(Tag.SOLVE, dt));
                            schwarzResid = 0.0;
                            convergesSchwarz = true;
                            for (int d = 0; d < nDoms; ++d) {
                                GmresResult results = (GmresResult) mpi.receive(d, Tag.GMRES_RESULT).getData();
                                if (results.converges == false) {
                                    convergesSchwarz = false;
                                    break;
                                }
                                schwarzResid += results.schwarzResid;
                            }
                            schwarzResid /= nDoms;

                            transferWatch.resume();
                            exchangeData(liteElems, mpi);
                            mpi.waitForAll(Tag.DATA_UPDATED);
                            transferWatch.suspend();

                            System.out.printf("     %d.   resid: %.2e\n", (schwarzIter + 1), schwarzResid);
                        }
                    }

                    if (!convergesNewton || !convergesSchwarz) {
                        mpi.sendAll(new MPIMessage(Tag.PREVIOUS_TIME_LEVEL));
                        state.cfl = cflObj.reduceCFL(state.cfl);
                        --state.steps;
                        LOG.warn("GMRES does not converge, CFL reduced to " + state.cfl);
                        continue outerloop;
                    }

                    mpi.sendAll(new MPIMessage(Tag.UPDATE_NEWTON));
                    mpi.waitForAll(Tag.NEWTON_UPDATED);
                }

                mpi.sendAll(new MPIMessage(Tag.LIMITER_and_NEXT_TIME_LEVEL, state.t + dt));
                double resid = 0.0;
                for (int d = 0; d < nDoms; ++d) {
                    resid += (double) mpi.receive(d, Tag.RESIDUUM).getData();
                }
                state.residuum = resid / nDoms;
                if (state.residuum == 0) {
                    LOG.error(" computation error ");
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
                state.transferTime = transferWatch.getTime();
                String info = infoToString(totalSteps, dt);
                if ((state.steps % par.saveRate) == 0) {
                    Solution solution = collectSolution(mpi);
                    if (par.animation) {
                        saveAnimationData(solution, state.steps);
                    }
                    saveData(solution);
                    LOG.info(info);
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
            state.transferTime = transferWatch.getTime();
            --state.steps;
            LOG.info("computation has finished in " + millisecsToTime(state.executionTime));
            LOG.info("primani: {}", tempWatch.toString());
            LOG.info("odesilani: {}", tempWatch2.toString());

            LOG.info("collecting data");
            Solution solution = collectSolution(mpi);

            return solution;
        } finally {
            mpi.sendAll(new MPIMessage(Tag.CLOSE));
            mpi.close();
        }
    }

    // vlastni vypocet
    public Solution localSolve() throws IOException {
        if (par.explicitTimeIntegration) {
            return localSolveExplicit();
        } else {
            return localSolveImplicit();
        }
    }

    public Solution localSolveImplicit() throws IOException {
        int nElems = elems.length;
        LinearSolver linSolver = LinearSolver.factory(par, elems, dofs);
        JacobiAssembler assembler = new JacobiAssembler(elems, par);
        double[] x = new double[dofs];
        StopWatch watch = new StopWatch();

        CFLSetup cflObj = new CFLSetup(par.cfl, par.varyCFL);
        double dto = -1;
        boolean converges = true;
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

            if (converges) {
                state.cfl = cflObj.getCFL(state.cfl, state.residuum);
            }
            // zastavovaci podminka
            if (state.cfl < (par.cfl / 20)) {
                LOG.error("algorithm does not converge - aborting computation");
                return mesh.getSolution();
            }

            // setting time step size - dt
            double dt = timeStep(state.cfl);
            if (dto == -1) {
                dto = dt;
            }
            if (state.t + dt > par.tEnd) {
                dt = par.tEnd - state.t;
            }

            //compute domain 
            if (par.solutionMonitorOn) {
                mesh.computeSolutionMonitor();
            }

            // compute artificial viscosity
            for (int i = 0; i < nElems; i++) {
                elems[i].limiter();
            }

            for (int s = 0; s < par.newtonIters; s++) {  // vnitrni iterace (Newton)
                // mesh deformation
                if (par.movingMesh) {
                    dfm.calculateForces(elems, dyn.getMeshMove());
                    dyn.computeBodyMove(dt, state.t, dfm.getFluidForces());
                    dfm.newMeshPosition(elems, par.orderInTime, dt, dto, dyn.getMeshMove());
                    if (state.t == 0) {
                        dfm.relaxFirstIteration(elems);
                    }
                    dfm.recalculateMesh(elems, par.order);
                }

                // matrix assembly
                assembler.assemble(dt, dto);

                // solve linear system
                Arrays.fill(x, 0.0);
                converges = linSolver.solve(x);

                if (!converges) {
                    getTimeLevelBack();
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
                    elems[i].updateW(x);
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
            // limiter
            for (int i = 0; i < nElems; i++) {
                elems[i].limitUnphysicalValues();
            }

            saveTimeLevel(state.t);

            mesh.updateTime(dt);
            if (par.movingMesh) {
                dfm.nextTimeLevel(elems);
                dyn.nextTimeLevel();
                dyn.savePositionsAndForces();
            }

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

    public Solution localSolveExplicit() throws IOException {
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
            ltsIter.iterate(state.t + dt);

            state.residuum = calculateResiduumW(dt);
            state.executionTime = watch.getTime();
            state.t += dt;
            saveResiduum(state.residuum, state.t, state.executionTime);
            saveTimeLevel(state.t);

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

    public void testDynamic(double dt) throws IOException {
        double t = 0;
        for (int step = 0; step <= par.steps; step++) {
            dyn.computeBodyMove(dt, t, new FluidForces(new double[2][dfm.nBodies], new double[1][dfm.nBodies], null, null, null));
            dyn.nextTimeLevel();
            dyn.savePositionsAndForces();
            t += dt;
            System.out.println(step + "-th iteration t = " + t);
        }
    }

    private void runFetcher(int nNodes, int fetcherPort, String masterIP, int masterPort) throws IOException {
        FetcherServer fetcher = new FetcherServer(fetcherPort, 3000, "matlab/pclist.txt");
        String args = "slave " + masterIP + " " + masterPort;
        ZipFile zip = new ZipFile("FlowPro.zip", "FlowPro.jar", args);
        fetcher.initFetcher(zip, nNodes);
        fetcher.start();
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
            if (par.order > 2) {
                Mat.save(sol.W, simulationPath + "animation/We" + (10000000 + step) + ".txt");
            }
            if (par.movingMesh) {
                Mat.save(sol.vertices, simulationPath + "animation/vertices" + (10000000 + step) + ".txt");
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
        private final Parameters par;
        double[] timeWCoef;
        double[] timeRHSCoef;

        public JacobiAssembler(Element[] elems, Parameters par) {
            this.elems = elems;
            this.par = par;
        }

        // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
        public void assemble(double dt, double dto) {  // , int newtonIter
            timeWCoef = tim.getWcoefficient();
            timeRHSCoef = tim.getRHScoefficient();
            switch (par.orderInTime) {
                case 1:
                    timeWCoef = new double[2];
                    timeWCoef[0] = 1 / dt;
                    timeWCoef[1] = -1 / dt;
                    break;
                case 2:
                    timeWCoef = new double[3];
                    timeWCoef[0] = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt);
                    timeWCoef[1] = -(dt + dto) / (dt * dto);  // -2/dt;
                    timeWCoef[2] = dt / (dto * (dt + dto));  // 1/(2*dt);
                default:
                    throw new RuntimeException("solver supports only first and second order in time");
            }

            AssemblerThread[] assemblers = new AssemblerThread[par.nThreads];

            // vlastni vypocet, parallelni beh
            for (int v = 0;
                    v < assemblers.length;
                    v++) {
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
                        elems[i].assembleJacobiMatrix(timeWCoef,timeRHSCoef);
                        elems[i].computeJacobiPreconditioner();
                    }
                }
            }
        }
    }

    public class LocalTimeStepIterator {

        public Element[] elems;
        private final int nEqs;
        private final Parameters par;

        public LocalTimeStepIterator(Element[] elems, Parameters par, double t) {
            this.elems = elems;
            this.par = par;
            nEqs = elems[0].getNEqs();
            for (Element elem : elems) {
                elem.tLTS = t;
                elem.tLTSold = t;
                int nBasis = elem.nBasis;
                elem.WLTS = new double[nBasis * nEqs];
                System.arraycopy(elem.W, 0, elem.WLTS, 0, nBasis * nEqs);
                switch (par.orderInTime) {
                    case 2:
                        elem.W1LTS = new double[nBasis * nEqs];
                        elem.W1LTSo = new double[nBasis * nEqs];
                        System.arraycopy(elem.W, 0, elem.W1LTS, 0, nBasis * nEqs);
                        System.arraycopy(elem.W, 0, elem.W1LTSo, 0, nBasis * nEqs);
                        break;
                    case 3:
                        elem.W1LTS = new double[nBasis * nEqs];
                        elem.W1LTSo = new double[nBasis * nEqs];
                        System.arraycopy(elem.W, 0, elem.W1LTS, 0, nBasis * nEqs);
                        System.arraycopy(elem.W, 0, elem.W1LTSo, 0, nBasis * nEqs);

                        elem.W2LTS = new double[nBasis * nEqs];
                        elem.W2LTSo = new double[nBasis * nEqs];
                        System.arraycopy(elem.W, 0, elem.W2LTS, 0, nBasis * nEqs);
                        System.arraycopy(elem.W, 0, elem.W2LTSo, 0, nBasis * nEqs);
                        break;
                }
            }
        }

        // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
        public void iterate(double tFinal) {  // , int newtonIter

            LocalTimeStepIteratorIteratorThread[] iteratorThrds = new LocalTimeStepIteratorIteratorThread[par.nThreads];

            // vlastni vypocet, parallelni beh
            for (int v = 0; v < iteratorThrds.length; v++) {
                iteratorThrds[v] = new LocalTimeStepIteratorIteratorThread(v, tFinal);
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

        private class LocalTimeStepIteratorIteratorThread extends Thread {

            private final int id;
            private final double tFinal;

            LocalTimeStepIteratorIteratorThread(int id, double tFinal) {
                this.id = id;
                this.tFinal = tFinal;
            }

            @Override
            public void run() {
                double CFLexp = par.cflLTS / par.order;
                boolean LTSdone = false;
                while (!LTSdone) {
                    LTSdone = true;
                    for (int i = id; i < elems.length; i += par.nThreads) {
                        if (elems[i].insideComputeDomain && elems[i].tLTS < tFinal) {
                            int s = 0;
                            for (int j = 0; j < elems[i].nFaces; j++) {
                                if (elems[i].TT[j] > -1) {
                                    if (elems[elems[i].TT[j]].tLTS >= elems[i].tLTS) {
                                        s++;
                                    }
                                } else {
                                    s++;
                                }
                            }
                            if (s == elems[i].nFaces) {
                                double dtLoc = elems[i].delta_t(CFLexp);
//                                if (elems[i].tLTS + dtLoc >= tFinal) {
//                                    dtLoc = tFinal - elems[i].tLTS;
//                                }
                                elems[i].computeExplicitStep(dtLoc);
                                elems[i].tLTSold = elems[i].tLTS;
                                elems[i].tLTS += dtLoc;
                            }
                        }

                        if (elems[i].insideComputeDomain && elems[i].tLTS < tFinal) {
                            LTSdone = false;
                        }
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
