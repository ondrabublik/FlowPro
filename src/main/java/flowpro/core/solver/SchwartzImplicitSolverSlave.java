package flowpro.core.solver;

import flowpro.core.*;
import flowpro.api.Mat;
import flowpro.core.parallel.*;
import flowpro.api.Equation;
import static flowpro.core.FlowProMain.*;
import flowpro.core.element.Element;
import flowpro.api.Dynamics;
import flowpro.api.MeshMove;
import flowpro.core.LinearSolvers.LinearSolver;
import flowpro.core.element.Implicit;
import flowpro.core.meshDeformation.*;
import litempi.*;

import java.io.*;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * tento program resi Navierovy-Stokesovy rovnice na nestrukturovane siti pomoci
 * nespojite Galerkinovy metody
 */
public class SchwartzImplicitSolverSlave extends SlaveSolver {

    private static final Logger LOG = LoggerFactory.getLogger(MasterSolver.class);

//    private final Object lock;
    // common parameters
    private final Deformation dfm;
    private final Dynamics dyn;
    private final Equation eqn;
    private final Parameters par;

    // master and loner parameters
    private Mesh[] meshes;
    private State state;
    private final Object lock;
    private String simulationPath;

    // slave and loner parameters
    private Mesh mesh;
    private Element[] elems;
    private final MPISlave mpiSlave;
    private int dofs;
//    private Solution solution;

    /**
     * Constructor for slave.
     *
     * @param slavePort
     * @throws IOException
     * @throws MPIException
     */
    public SchwartzImplicitSolverSlave(int slavePort) throws IOException, MPIException {
        mpiSlave = new MPISlave(slavePort);
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

    private void computeArtificialViscosity(boolean isFirstIter) {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.limiter(isFirstIter);
            }
        }
    }

    private void copyWo2W() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.copyWo2W();
            }
        }
    }

    private void copyW2Wo() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.copyW2Wo();
            }
        }
    }

    private void updateRHS(double x[]) {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                ((Implicit) elem.ti).updateRHS(x);
            }
        }
    }

    private void updateW(double x[]) {
        for (Element elem : elems) {
            elem.updateW(x);
        }
    }

    private void nextTimeLevelMassMatrixes() {
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.nextTimeLevelMassMatrixes();
            }
        }
    }

    private double calculateResiduumW(double dt) {
        double resid = 0;
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                resid += elem.calculateResiduumW(dt);
            }
        }

        return resid / elems.length;
    }

    public void solve() throws IOException, MPIException {
        LOG.info("slave is running...");
        MPISlave mpi = mpiSlave;
        mesh = meshes[0];
        try {
            LinearSolver linSolver = LinearSolver.factory(elems, par);
            JacobiAssembler assembler = new JacobiAssembler(elems, par);
            double[] x = new double[dofs];
            double[] y = new double[dofs];
            double dt = -1.0;
            double dto;
            boolean isFirstIter = true;

            while (true) {
                /* receive message */
                MPIMessage inMsg = mpi.receive();
//                LOG.debug("received {}", inMsg.tag);
                MPIMessage outMsg = null;

                switch (inMsg.tag) {
                    case Tag.TIME_STEP_REQ: // vypocet dt
                        double CFL = (double) inMsg.getData();
                        double locDt = timeStep(CFL);
                        computeArtificialViscosity(isFirstIter);

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
                        dfm.newMeshPositionAndVelocity(elems, elems[0].ti.getOrder(), disp.getDto(), disp.getDt(), disp.getMeshMove());
                        if (isFirstIter) {
                            dfm.relaxFirstIteration(elems, disp.getDt());
                        }
                        dfm.recalculateMesh(elems);
                        outMsg = new MPIMessage(Tag.MESH_POSITION_UPDATED);
                        break;

                    case Tag.ALE_NEXT_TIME_LEVEL:
                        dfm.nextTimeLevel(elems);
                        nextTimeLevelMassMatrixes();
                        outMsg = new MPIMessage(Tag.NEW_MESH_POSITION_UPDATED);
                        break;

                    case Tag.ASSEMBLE_AND_SOLVE: // inicialization if equation Ax=b
                        // dalo by se tomu vyhnout kdybysme ukladali souboru dt a dto
                        if (isFirstIter) {  //  if dt has not yet been set
                            dt = (double) inMsg.getData();
                            dto = dt;
                        } else {  // if dt has already been set
                            dto = dt;
                            dt = (double) inMsg.getData();
                        }
                        // set equation object state
                        eqn.setState(mesh.t + dt, dt);

                        assembler.assemble(dt, dto);
                        Arrays.fill(y, 0.0);
                    // NO BREAK! continue to the following tag                    // NO BREAK! continue to the following tag

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
                        double residuum = calculateResiduumW(dt);
                        copyW2Wo();
                        mesh.updateTime(dt);
                        outMsg = new MPIMessage(Tag.RESIDUUM, residuum);
                        isFirstIter = false;
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
                                    yElem[ind] = y[elems[j].gIndex[ind]];
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
                                    y[elems[j].gIndex[ind]] = dataReceive1.y[ind];
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
                        copyWo2W();
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

    public void saveResiduum(double residuum, double t, double CPU) {
        try {
            FileWriter fw;
            fw = new FileWriter(simulationPath + "residuum.txt", true);
            BufferedWriter bw = new BufferedWriter(fw);
            try (PrintWriter out = new PrintWriter(bw)) {
                out.println(Double.toString(residuum) + " " + Double.toString(t) + " " + Double.toString(CPU));
            }
        } catch (IOException ex) {
            LOG.warn("error occured while writing to a file: " + ex.getMessage());
        }
    }
}
