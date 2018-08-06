//package flowpro.core;
//
//import flowpro.core.meshDeformation.*;
//import flowpro.api.FlowProProperties;
//import flowpro.api.Mat;
//import flowpro.core.parallel.GmresResult;
//import flowpro.core.parallel.Domain;
//import flowpro.core.parallel.Tag;
//import flowpro.core.parallel.LiteElement;
//import flowpro.core.iterativeLinearSolvers.IterativeLinearSolver;
//import flowpro.user.dynamics.Rigid2D;
//import flowpro.api.Equation;
//import static flowpro.core.FlowProMain.*;
//import flowpro.core.parallel.Domain.Subdomain;
//import flowpro.core.Mesh.Element;
//import fetcher.FetcherServer;
//import fetcher.ZipFile;
//import flowpro.api.FluidForces;
//import litempi.*;
//
//import java.io.*;
//import java.util.Arrays;
//import java.util.logging.Level;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;
//import org.apache.commons.lang3.time.StopWatch;
//
///**
// * tento program resi Navierovy-Stokesovy rovnice na nestrukturovane siti pomoci
// * nespojite Galerkinovy metody
// */
//public class SolverExperimental {
//
//    private static final Logger LOG = LoggerFactory.getLogger(Solver.class);
//    private static final StopWatch tempWatch = new StopWatch();
//    private static final StopWatch tempWatch2 = new StopWatch();
//
////    private final Object lock;
//    // common parameters
//    private final Deformation dfm;
//    private final Rigid2D dyn;
//    private final Equation eqn;
//    private final Parameters par;
//
//    // master and loner parameters
//    private Mesh[] meshes;
//    private State state;
//    private Domain domain;
//    private final Object lock;
//    private String simulationPath;
//
//    // slave and loner parameters
//    private Mesh mesh;
//    private Element[] elems;
//    private MPISlave mpiSlave;
//    private int dofs;
////    private Solution solution;
//
//    /**
//     * Constructor for master and loner.
//     *
//     * @param simulationPath
//     * @param meshes
//     * @param dyn
//     * @param eqn
//     * @param par
//     * @param state
//     * @param domain
//     * @param lock
//     */
//    public SolverExperimental(String simulationPath, Mesh[] meshes, Rigid2D dyn,
//            Equation eqn, Parameters par, State state, Domain domain, Object lock) {
//        this.simulationPath = simulationPath;
//        this.meshes = meshes;
//        this.dyn = dyn;
//        this.eqn = eqn;
//        this.par = par;
//        this.state = state;
//        this.domain = domain;
//        mesh = meshes[0];
//        elems = mesh.getElems();
//        this.dofs = mesh.dofs;
//        this.dfm = mesh.getDfm();
//        this.lock = lock;
//    }
//
//    /**
//     * Constructor for slave.
//     *
//     * @param masterIP
//     * @param masterPort
//     * @throws IOException
//     * @throws MPIException
//     */
//    public SolverExperimental(String masterIP, int masterPort) throws IOException, MPIException {
//        mpiSlave = new MPISlave(masterIP, masterPort);
//        MPIMessage msg = mpiSlave.receive();
//
//        if (msg.tag != Tag.INIT_DATA) {
//            throw new MPIException("expected InitMessage");
//        }
//
//        mesh = (Mesh) msg.getData();
//        meshes = new Mesh[]{mesh};
//        mesh.init();
//        elems = mesh.getElems();
//
//        par = mesh.getPar();
//        eqn = mesh.getEqn();
//        dfm = mesh.getDfm();
//        dyn = null;
//        dofs = mesh.dofs;
//        lock = null;
//
//        mpiSlave.send(new MPIMessage(Tag.DATA_INITIALIZED));
//    }
//
//    private double timeStep(double CFL) {
//        double dt = Double.MAX_VALUE;
//        for (Element elem : elems) {
//            if (elem.insideComputeDomain) {
//                double loc_dt = elem.delta_t(CFL);
//                if (loc_dt < dt) {
//                    dt = loc_dt;
//                }
//            }
//        }
//
//        return dt;
//    }
//
//    private void callLimiters(String set) {
//        switch (set) {
//            case ("all"):
//                for (Element elem : elems) {
//                    if (elem.insideComputeDomain) {
//                        elem.limitUnphysicalValues();
//                        elem.limiter();
//                    }
//                }
//            case ("interior"):
//                for (Element elem : elems) {
//                    if (elem.insideComputeDomain && elem.insideMetisDomain) {
//                        elem.limitUnphysicalValues();
//                        elem.limiter();
//                    }
//                }
//            case ("overlap"):
//                for (Element elem : elems) {
//                    if (elem.insideComputeDomain && !elem.insideMetisDomain) {
//                        elem.limitUnphysicalValues();
//                        elem.limiter();
//                    }
//                }
//                break;
//        }
//    }
//
//    private void copyWo2W() {
//        for (Element elem : elems) {
//            if (elem.insideComputeDomain) {
//                elem.copyWo2W();
//            }
//        }
//    }
//
//    private void copyW2Wo(String set) {
//        switch (set) {
//            case ("all"):
//                for (Element elem : elems) {
//                    if (elem.insideComputeDomain) {
//                        elem.copyW2Wo();
//                    }
//                }
//                break;
//            case ("interior"):
//                for (Element elem : elems) {
//                    if (elem.insideComputeDomain && elem.insideMetisDomain) {
//                        elem.copyW2Wo();
//                    }
//                }
//                break;
//            case ("overlap"):
//                for (Element elem : elems) {
//                    if (elem.insideComputeDomain && !elem.insideMetisDomain) {
//                        elem.copyW2Wo();
//                    }
//                }
//                break;
//        }
//    }
//
//    /**
//     * b = b - Ax
//     *
//     * @param x
//     */
//    private void updateRHS(double x[]) {
//        for (Element elem : elems) {
//            if (elem.insideComputeDomain) {
//                elem.updateRHS(x);
//            }
//        }
//    }
//
//    private void updateW(double x[], String set) {
//        switch (set) {
//            case ("all"):
//                for (Element elem : elems) {
//                    elem.updateW(x);
//                }
//                break;
//            case ("interior"):
//                for (Element elem : elems) {
//                    if (elem.insideMetisDomain) {
//                        elem.updateW(x);
//                    }
//                }
//                break;
//            case ("overlap"):
//                for (Element elem : elems) {
//                    if (!elem.insideMetisDomain) {
//                        elem.updateW(x);
//                    }
//                }
//                break;
//        }
//    }
//
//    /**
//     *
//     * @param dt
//     * @return L1norm(W - Wo)
//     */
//    private double calculateResiduumW(double dt) {
//        double resid = 0;
//        for (Element elem : elems) {
//            if (elem.insideMetisDomain) {
//                resid += elem.calculateResiduumW(dt);
//            }
//        }
//        return resid / elems.length;
//    }
//
//    public void slaveSolve() throws IOException, MPIException, InterruptedException {
//        LOG.info("slave is running...");
//        MPISlave mpi = mpiSlave;
//        mesh = meshes[0];
//
//        SlaveSolverParameters slavePar = new SlaveSolverParameters();
//        slavePar.linSolver = IterativeLinearSolver.factory(par, elems, dofs);
//        slavePar.assembler = new JacobiAssembler(elems, par);
//        slavePar.x = new double[dofs];
//        slavePar.y = new double[dofs];
//        slavePar.dt = -1.0;
//        slavePar.dto = -1.0;
//        slavePar.firstTimeStep = true;
//        slavePar.mpi = mpi;
//
//        while (true) {
//            /* receive message */
//            MPIMessage inMsg = mpi.receive();
//            if(inMsg.tag == Tag.CLOSE){
//                LOG.info("received command to shut down");
//                break;
//            }
//            SlaveSolveThread sst = new SlaveSolveThread(slavePar, inMsg);
//            sst.start();
//        }
//    }
//
//    public class SlaveSolverParameters {
//
//        IterativeLinearSolver linSolver;
//        JacobiAssembler assembler;
//        double[] x;
//        double[] y;
//        double dt;
//        double dto;
//        boolean firstTimeStep;
//        MPISlave mpi;
//    }
//
//    public class SlaveSolveThread extends Thread {
//
//        SlaveSolverParameters sp;
//        MPIMessage inMsg;
//
//        SlaveSolveThread(SlaveSolverParameters sp, MPIMessage inMsg) {
//            this.sp = sp;
//            this.inMsg = inMsg;
//        }
//
//        @Override
//        public void run() {
//            MPIMessage outMsg = null;
//            switch (inMsg.tag) {
//                case Tag.TIME_STEP_REQ: // vypocet dt
//                    double CFL = (double) inMsg.getData();
//                    double locDt = timeStep(CFL);
//                    outMsg = new MPIMessage(Tag.TIME_STEP, locDt);
//                    break;
//
//                case Tag.LIMITERS:
//                    String set = (String) inMsg.getData();
//                    callLimiters(set);
//                    break;
//
//                case Tag.ALE_CALCULATE_FORCES:
//                    dfm.calculateForces(elems);
//                    outMsg = new MPIMessage(Tag.FORCES, new ForcesAndDisplacements(dfm.getFluidForces()));
//                    break;
//
//                case Tag.ALE_NEW_MESH_POSITION:
//                    ForcesAndDisplacements disp = (ForcesAndDisplacements) inMsg.getData();
//                    dfm.newMeshPosition(elems, par.orderInTime, disp.getDto(), disp.getDt(), disp.getMeshMove());
//                    dfm.recalculateMesh(elems, par.order);
//                    outMsg = new MPIMessage(Tag.MESH_POSITION_UPDATED);
//                    break;
//
//                case Tag.ASSEMBLE_INTERIOR: // interior inicialization if equation Ax=b
//                    // dalo by se tomu vyhnout kdybysme ukladali souboru dt a dto
//                    if (sp.firstTimeStep) {  //  if dt has not yet been set
//                        sp.dt = (double) inMsg.getData();
//                        sp.dto = sp.dt;
//                        sp.firstTimeStep = false;
//                    } else {  // if dt has already been set
//                        sp.dto = sp.dt;
//                        sp.dt = (double) inMsg.getData();
//                    }
//                    sp.assembler.assemble(sp.dt, sp.dto, "interior");
//                    outMsg = new MPIMessage(Tag.ASSEMBELED);
//                    break;
//
//                case Tag.ASSEMBLE_OVERLAP_and_SOLVE: // interior inicialization if equation Ax=b
//                    // assembling overlaping part of Jacobi matrix
//                    sp.assembler.assemble(sp.dt, sp.dto, "overlap");
//                    Arrays.fill(sp.y, 0.0);
//                // no break !!!!!!!!!!!!!
//
//                case Tag.SOLVE: // solving linear equation Ax = b - Ay
//                    Arrays.fill(sp.x, 0.0);
//                    boolean converges = sp.linSolver.solve(sp.x);
//                    Mat.plusEqual(sp.y, sp.x);
//                    updateRHS(sp.x);  // b = b - Ax
//                    double schwarzResid = Mat.L1Norm(sp.x) / sp.x.length;
//
//                    outMsg = new MPIMessage(Tag.GMRES_RESULT, new GmresResult(converges, schwarzResid));
//                    break;
//
//                case Tag.INTERIOR_UPDATE: // update newton
//                    updateW(sp.y, "interior");
//                    outMsg = new MPIMessage(Tag.INTERIOR_UPDATED);
//                    break;
//
//                case Tag.INTERIOR_LIMITER_and_NEXT_TIME_LEVEL: // next time level               
//                    double residuum = calculateResiduumW(sp.dt);
//                    callLimiters("interior");
//                    copyW2Wo("interior");
//                    mesh.updateTime(sp.dt);
//                    outMsg = new MPIMessage(Tag.RESIDUUM, residuum);
//                    break;
//
//                case Tag.OVERLAP_UPDATE_and_LIMITER_and_NEXT_TIME: // interior inicialization if equation Ax=b
//                    // updating newton iteration in overlapping part
//                    updateW(sp.y, "overlap");
//
//                    // next time level for overlaping part
//                    int newtonIter = (int) inMsg.getData();
//                    if (newtonIter == 0) {
//                        callLimiters("overlap");
//                        copyW2Wo("overlap");
//                    }
//                    outMsg = new MPIMessage(Tag.OVERLAP_UPDATED);
//                    break;
//
//                case Tag.DATA_REQ: // saving data into central structure
//                    int nSave = mesh.save.length;
//                    LiteElement[] dataSend = new LiteElement[nSave];
//                    for (int i = 0; i < nSave; i++) {
//                        int j = mesh.save[i];
//                        double[] yElem = new double[mesh.nEqs * elems[j].nBasis];
//                        for (int m = 0; m < mesh.nEqs; m++) {
//                            for (int p = 0; p < elems[j].nBasis; p++) {
//                                int ind = elems[j].nBasis * m + p;
//                                yElem[ind] = sp.y[elems[j].gi_U[ind]];
//                            }
//                        }
//                        dataSend[i] = new LiteElement(j, yElem);
//                    }
//                    outMsg = new MPIMessage(Tag.DATA_SLAVE_TO_MASTER, dataSend);
//                    break;
//
//                case Tag.DATA_MASTER_TO_SLAVE: // loading data from central structure
//                    LiteElement[] dataReceive = (LiteElement[]) inMsg.getData();
//                    for (LiteElement dataReceive1 : dataReceive) {
//                        int j = dataReceive1.index;
//                        for (int m = 0; m < mesh.nEqs; m++) {
//                            for (int p = 0; p < elems[j].nBasis; p++) {
//                                int ind = elems[j].nBasis * m + p;
//                                sp.y[elems[j].gi_U[ind]] = dataReceive1.y[ind];
//                            }
//                        }
//                    }
//                    outMsg = new MPIMessage(Tag.DATA_UPDATED);
//                    break;
//
//                case Tag.SOLUTION_REQ: // zaverecne odeslani vsech dat (krome vnejsiho prekrivu)
//                    Solution localSolution = mesh.getSolution();
////                        Solution localSolution = new Solution(mesh.getW(), mesh.getAvgW(), mesh.getDetailW());
//                    LOG.info("sending solution");
//                    outMsg = new MPIMessage(Tag.SOLUTION, localSolution);
//                    break;
//
//                case Tag.PREVIOUS_TIME_LEVEL: // GMRES does not converge, go to previous time level
//                    copyWo2W();
//                    break;
//
//                case Tag.RESET_OUTPUT_WRITER:  // reset output stream
//                    try {
//                        sp.mpi.reset();
//                    } catch (MPIException ex) {
//                        java.util.logging.Logger.getLogger(Solver.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                    break;
//
//                case Tag.GET_MEMORY:
//                    long total = Runtime.getRuntime().totalMemory();
//                    long free = Runtime.getRuntime().freeMemory();
//                    long max = Runtime.getRuntime().maxMemory();
//                    long usage = total - free;
//                    String output = "Usage = " + usage / 1048576 + "Mb (" + (int) ((double) usage / max * 100) + ")%";
//                    outMsg = new MPIMessage(Tag.MEMORY_USAGE, output);
//                    break;
//            }
//
//            /* send message */
//            if (outMsg != null) {
//                try {
//                    sp.mpi.send(outMsg);
//                } catch (MPIException ex) {
//                    java.util.logging.Logger.getLogger(Solver.class.getName()).log(Level.SEVERE, null, ex);
//                }
//            }
//        }
//    }
//
//    private String infoToString(int totalSteps, double dt) throws IOException {
//        String timeStr = millisecsToTime(state.getOverallExecutionTime());
//
//        return String.format("%d/%d  resid: %.2e,  dt: %.1e,  t: %.2f,  CFL: %1.2f,  CPU: %s",
//                state.steps, totalSteps, state.residuum, dt, state.t, state.cfl, timeStr);
//    }
//
//    private Solution collectSolution(MPIMaster mpi) throws MPIException {
//        mpi.sendAll(new MPIMessage(Tag.SOLUTION_REQ));
//        Solution[] sols = new Solution[mpi.nSlaves];
//        for (int d = 0; d < mpi.nSlaves; ++d) {
//            sols[d] = (Solution) mpi.receive(d, Tag.SOLUTION).getData();
//        }
//        return new Solution(sols, domain);
//    }
//
//    private void exchangeData(LiteElement[] liteElems, MPIMaster mpi) throws MPIException {
//        // downloading data from the slave nodes into the central structure
//        mpi.sendAll(new MPIMessage(Tag.DATA_REQ));
//        for (int d = 0; d < domain.nDoms; ++d) {
//            int[] mapL2G = domain.getSubdomain(d).mapL2G;
//            LiteElement[] dataRcv = (LiteElement[]) mpi.receive(d, Tag.DATA_SLAVE_TO_MASTER).getData();
//            for (LiteElement dataRcv1 : dataRcv) {
//                //liteElems[mapL2G[dataRcv1.index]] = dataRcv1;
//                liteElems[mapL2G[dataRcv1.index]] = new LiteElement(dataRcv1.index, dataRcv1.y);
//            }
//        }
//
//        // uploading data from the central structure onto the slave nodes
//        for (int d = 0; d < domain.nDoms; ++d) {
//            Subdomain subdom = domain.getSubdomain(d);
//            int[] mapG2L = subdom.mapG2L;
//            int[] load = subdom.load;
//            LiteElement[] dataSend = new LiteElement[load.length];
//            for (int i = 0; i < load.length; i++) {
//                dataSend[i] = new LiteElement(mapG2L[load[i]], liteElems[load[i]].y);
//            }
//            mpi.send(new MPIMessage(Tag.DATA_MASTER_TO_SLAVE, dataSend), d);
//        }
//    }
//
//    public Solution masterSolve() throws MPIException, IOException {
//        int nDoms = domain.nDoms;
//        runFetcher(nDoms, par.fetcherPort, par.masterIP, par.masterPort);
//        MPIMaster mpi = new MPIMaster(nDoms, par.masterIP, par.masterPort);
//        StopWatch watch = new StopWatch();
//        try {
//            LOG.info("sending initial data");
//            for (int d = 0; d < nDoms; ++d) {
//                mpi.send(new MPIMessage(Tag.INIT_DATA, meshes[d]), d);
//                meshes[d] = null;  // garbage collector should clear the memory
//            }
//            meshes = null;  // garbage collector should clear the memory
//            mpi.waitForAll(Tag.DATA_INITIALIZED);
//
//            // central structure
//            LiteElement[] liteElems = new LiteElement[domain.nElems];
//            boolean converges = true;
//            int totalSteps = state.steps + par.steps;
//
//            LOG.info("computation has started...");
//            watch.start();
//            outerloop:
//            // main loop ============================== 
//            for (++state.steps; state.steps <= totalSteps && state.residuum > par.residuum
//                    && state.t < par.tEnd; ++state.steps) {
//                if (converges) {
//                    state.cfl += par.cfl / 20;
//                    if (state.cfl > par.cfl) {
//                        state.cfl = par.cfl;
//                    }
//                }
//                // stop condition
//                if (state.cfl < (par.cfl / 20)) {
//                    LOG.error("algorithm does not converge - aborting computation");
//                    System.exit(1);
//                }
//
//                // time step computation
//                mpi.sendAll(new MPIMessage(Tag.TIME_STEP_REQ, state.cfl));
//                double dt = Mat.min(mpi.receiveAllDouble(Tag.TIME_STEP));
//
//                // mesh deformation
//                if (par.movingMesh) {
//                    mpi.sendAll(new MPIMessage(Tag.ALE_CALCULATE_FORCES));
//                    ForcesAndDisplacements[] forDis = new ForcesAndDisplacements[nDoms];
//                    for (int d = 0; d < nDoms; ++d) {
//                        forDis[d] = (ForcesAndDisplacements) mpi.receive(d, Tag.FORCES).getData();
//                    }                
//                    dyn.computeBodyMove(dt, state.t, forDis[0].combine(forDis));
//                    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//                    mpi.sendAll(new MPIMessage(Tag.ALE_NEW_MESH_POSITION, new ForcesAndDisplacements(dt, -100000000000000000.0 ,dyn.getMeshMove()))); 
//                    mpi.waitForAll(Tag.MESH_POSITION_UPDATED);
//                }
//
//                // newton loop ============================== 
//                for (int newtonIter = 0; newtonIter < par.newtonIters; newtonIter++) {
//                    mpi.sendAll(new MPIMessage(Tag.ASSEMBLE_INTERIOR, dt));
//                    exchangeData(liteElems, mpi);
//                    mpi.waitForAll(Tag.DATA_UPDATED);
//                    mpi.sendAll(new MPIMessage(Tag.OVERLAP_UPDATE_and_LIMITER_and_NEXT_TIME, newtonIter));
//                    mpi.waitForAll(Tag.OVERLAP_UPDATED);
//                    mpi.waitForAll(Tag.ASSEMBELED); // wait for assembler
//                    mpi.sendAll(new MPIMessage(Tag.ASSEMBLE_OVERLAP_and_SOLVE));
//
//                    converges = true;
//                    for (int d = 0; d < nDoms; ++d) {
//                        GmresResult results = (GmresResult) mpi.receive(d, Tag.GMRES_RESULT).getData();
//                        if (results.converges == false) {
//                            converges = false;
//                        }
//                    }
//
//                    // schwartz loop ============================== 
//                    double schwarzResid = Double.MAX_VALUE;
//                    for (int schwarzIter = 0; schwarzIter < par.schwarzIters && schwarzResid > par.schwarzTol; ++schwarzIter) {
//
//                        exchangeData(liteElems, mpi);
//                        mpi.waitForAll(Tag.DATA_UPDATED);
//
//                        mpi.sendAll(new MPIMessage(Tag.SOLVE, dt));
//                        schwarzResid = 0.0;
//                        converges = true;
//                        for (int d = 0; d < nDoms; ++d) {
//                            GmresResult results = (GmresResult) mpi.receive(d, Tag.GMRES_RESULT).getData();
//                            if (results.converges == false) {
//                                converges = false;
//                            }
//                            schwarzResid += results.schwarzResid;
//                        }
//                        schwarzResid /= nDoms;
//
//                        if (!converges) {
//                            mpi.sendAll(new MPIMessage(Tag.PREVIOUS_TIME_LEVEL));
//                            state.cfl /= 2;
//                            --state.steps;
//                            LOG.warn("GMRES does not converge, CFL reduced to " + state.cfl);
//                            continue outerloop;
//                        }
//
//                        System.out.printf("     %d.   resid: %.2e\n", schwarzIter + 1, schwarzResid);
//                    }
//
//                    mpi.sendAll(new MPIMessage(Tag.INTERIOR_UPDATE));
//                    mpi.waitForAll(Tag.INTERIOR_UPDATED);
//                }
//
//                // next time level
//                mpi.sendAll(new MPIMessage(Tag.INTERIOR_LIMITER_and_NEXT_TIME_LEVEL));
//                state.residuum = Mat.sum(mpi.receiveAllDouble(Tag.RESIDUUM)) / nDoms;
//
//                if (par.movingMesh) {
//                    dyn.savePositionsAndForces();
//                }
//
//                state.t += dt;
//                state.executionTime = watch.getTime();
//                String info = infoToString(totalSteps, dt);
//                if ((state.steps % par.saveRate) == 0 && state.steps > 1) {
//                    Solution solution = collectSolution(mpi);
//                    if (par.movingMesh) {
//                        saveAnimationData(solution, state.steps);
//                    } else {
//                        saveData(solution);
//                    }
//                    LOG.info(info);
//                }
//                System.out.printf(info);
//                System.out.println();
//
//                // free memory
//                if (state.steps % 20 == 0) {
//                    mpi.reset();
//                    mpi.sendAll(new MPIMessage(Tag.RESET_OUTPUT_WRITER));
//                }
//            }
//            watch.stop();
//            state.executionTime = watch.getTime();
//            state.transferTime = 0;
//            --state.steps;
//            LOG.info("computation has finished in " + millisecsToTime(state.executionTime));
//            LOG.info("primani: {}", tempWatch.toString());
//            LOG.info("odesilani: {}", tempWatch2.toString());
//
//            LOG.info("collecting data");
//            Solution solution = collectSolution(mpi);
//
//            return solution;
//        } finally {
//            mpi.sendAll(new MPIMessage(Tag.CLOSE));
//            mpi.close();
//        }
//    }
//
//    private void runFetcher(int nNodes, int fetcherPort, String masterIP, int masterPort) throws IOException {
//        FetcherServer fetcher = new FetcherServer(fetcherPort, 3000, "matlab/pclist.txt");
//        String args = "slave " + masterIP + " " + masterPort;
//        ZipFile zip = new ZipFile("DGFEM2D.zip", "DGFEM2D.jar", args);
//        fetcher.initFetcher(zip, nNodes);
//        fetcher.start();
//    }
//
//    public void saveData(Solution sol) throws IOException {
//        synchronized (lock) {
//            state.save();
//            eqn.saveReferenceValues(simulationPath + REF_VALUE_FILE_NAME);
//            Mat.save(sol.avgW, simulationPath + "W.txt");
//            Mat.save(sol.W, simulationPath + "We.txt");
//            // Mat.save(sol.detailW, simulationPath + "Wdetail.txt");
//            if (par.movingMesh) {
//                Mat.save(sol.vertices, simulationPath + "PXY.txt");
//            }
//
////            Mat.save(sol.detailW, simulationPath + "Wdetail.txt");
////            Mat.save(mesh.getViscosity(), simulationPath + "artificialViscosity.txt");
////            sol.saveSolution(simulationPath + "results.txt");
//            lock.notify();
//        }
//        LOG.info("results have been saved into " + simulationPath);
//    }
//
//    public void saveAnimationData(Solution sol, int step) throws IOException {
//        File directory = new File(simulationPath + "animation");
//        if (!directory.exists()) {
//            directory.mkdir();
//        }
//        synchronized (lock) {
//            Mat.save(sol.avgW, simulationPath + "animation/W" + (10000000 + step) + ".txt");
//            Mat.save(sol.vertices, simulationPath + "animation/PXY" + (10000000 + step) + ".txt");
//            lock.notify();
//        }
//        LOG.info("results have been saved into " + simulationPath);
//    }
//
//    public static class State {
//
//        public final String stateFilePath;
//        public final int order; // order of accuracy in space
//        private boolean hasOrderChanged;
//
//        int steps; // number of time steps that has already been taken
//        double t;  // physical time
//        double cfl;
//        double residuum;
//        long executionTime;  // current execution time
//        long transferTime;   // current transfer time
//
//        long initExecutionTime; // initial execution time
//        long initTransferTime;  // initial transfer time
//
//        public State(String stateFilePath, int order) throws IOException {
//            this.stateFilePath = stateFilePath;
//            this.order = order;
//            hasOrderChanged = false;
//
//            steps = 0;
//            t = 0.0;
//            cfl = 0.0;
//            residuum = Double.MAX_VALUE;
//            initExecutionTime = 0;
//            executionTime = 0;
//            initExecutionTime = 0;
//            transferTime = 0;
//        }
//
//        public void load() throws IOException {
//            try {
//                FlowProProperties stateProperties = new FlowProProperties();
//                stateProperties.load(new FileInputStream(stateFilePath));
//
//                t = stateProperties.getDouble("t");
//                steps = stateProperties.getInt("steps");
//                initExecutionTime = stateProperties.getLong("CPU");
//                initTransferTime = stateProperties.getLong("transfer");
//                cfl = stateProperties.getDouble("CFL");
//
//                int oldOrder = stateProperties.getInt("order");
//                if (oldOrder != order) {
//                    hasOrderChanged = true;
//                }
//            } catch (IOException ex) {
//                throw new IOException("file " + stateFilePath + " has a wrong format: "
//                        + ex.getMessage());
//            }
//        }
//
//        public void save() throws IOException {
//            FlowProProperties output = new FlowProProperties();
//
//            output.setProperty("t", Double.toString(t));
//            output.setProperty("steps", Integer.toString(steps));
//            output.setProperty("CPU", Long.toString(getOverallExecutionTime()));
//            output.setProperty("transfer", Long.toString(getOverallTransferTime()));
//            output.setProperty("CFL", Double.toString(cfl));
//            output.setProperty("residuum", Double.toString(residuum));
//
//            output.setProperty("order", Integer.toString(order));
//
//            output.store(new FileOutputStream(stateFilePath), null);
//        }
//
//        public boolean hasOrderChanged() {
//            return hasOrderChanged;
//        }
//
//        public long getOverallExecutionTime() {
//            return initExecutionTime + executionTime;
//        }
//
//        public long getOverallTransferTime() {
//            return initTransferTime + transferTime;
//        }
//    }
//
//    public class JacobiAssembler {
//
//        public Element[] elems;
//        private int nEqs;
//        private final Parameters par;
//        private double[] a1;
//        private double[] a2;
//        private double[] a3;
//        private double[] dual;
//        private double[] coeffsPhys;
//        private double[] coeffsDual;
//
//        public JacobiAssembler(Element[] elems, Parameters par) {
//            this.elems = elems;
//            this.par = par;
//            nEqs = elems[0].getNEqs();
//        }
//
//        // vytvoreni vlaken, paralelni sestaveni lokalnich matic a plneni globalni matice
//        public void assemble(double dt, double dto, String set) {  // , int newtonIter
//            if ("singleTime".equals(par.timeMethod)) {
//                coeffsPhys = new double[nEqs]; // user defined
//                coeffsDual = new double[nEqs];
//                for (int i = 0; i < nEqs; i++) {
//                    coeffsPhys[i] = 1;
//                    coeffsDual[i] = 0;
//                }
//            } else {
//                coeffsPhys = par.coeffsPhys;
//                coeffsDual = par.coeffsDual;
//            }
//
//            a1 = new double[nEqs];
//            a2 = new double[nEqs];
//            a3 = new double[nEqs];
//            dual = new double[nEqs];
//            if (par.orderInTime == 1) {
//                for (int i = 0; i < nEqs; i++) {
//                    a1[i] = coeffsPhys[i] / dt;
//                    a2[i] = -coeffsPhys[i] / dt;
//                    a3[i] = 0.0;
//                }
//            } else if (par.orderInTime == 2) {
//                for (int i = 0; i < nEqs; i++) {
//                    a1[i] = coeffsPhys[i] * (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
//                    a2[i] = -coeffsPhys[i] * (dt + dto) / (dt * dto);  // -2/dt;
//                    a3[i] = coeffsPhys[i] * dt / (dto * (dt + dto));  // 1/(2*dt);
//                }
//            } else {
//                throw new RuntimeException("solver supports only first and second order in time");
//            }
//            for (int i = 0; i < nEqs; i++) {
//                dual[i] = coeffsDual[i] / dt;
//            }
//
//            AssemblerThread[] assemblers = new AssemblerThread[par.nThreads];
//
//            // vlastni vypocet, parallelni beh
//            for (int v = 0; v < assemblers.length; v++) {
//                assemblers[v] = new AssemblerThread(v, set);
//                assemblers[v].start();
//            }
//
//            try {
//                for (AssemblerThread assembler : assemblers) {
//                    assembler.join();
//                }
//            } catch (java.lang.InterruptedException e) {
//                System.err.println(e);
//                System.exit(1);
//            }
//        }
//
//        private class AssemblerThread extends Thread {
//
//            private final int id;
//            private boolean overlapAssemblerSet;
//            private boolean interiorAssemblerSet;
//
//            AssemblerThread(int id, String set) {
//                this.id = id;
//                switch (set) {
//                    case "all":
//                        overlapAssemblerSet = true;
//                        interiorAssemblerSet = true;
//                        break;
//                    case "overlap":
//                        overlapAssemblerSet = true;
//                        interiorAssemblerSet = false;
//                        break;
//                    case "interior":
//                        overlapAssemblerSet = false;
//                        interiorAssemblerSet = true;
//                        break;
//                }
//            }
//
//            @Override
//            public void run() {
//                for (int i = id; i < elems.length; i += par.nThreads) {
//                    if (interiorAssemblerSet && elems[i].insideComputeDomain && elems[i].insideAssemblerDomain) {
//                        elems[i].assembleJacobiMatrix(a1, a2, a3, dual);
//                        elems[i].computeJacobiPreconditioner();
//                    }
//                    if (overlapAssemblerSet && elems[i].insideComputeDomain && !elems[i].insideAssemblerDomain) {
//                        elems[i].assembleJacobiMatrix(a1, a2, a3, dual);
//                        elems[i].computeJacobiPreconditioner();
//                    }
//                }
//            }
//        }
//    }
//}
