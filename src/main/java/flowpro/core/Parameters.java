package flowpro.core;

import flowpro.api.DomainTransformationObject;
import flowpro.api.FlowProProperties;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.Serializable;
import java.net.URL;

/**
 *
 * @author ales
 */
public class Parameters implements Serializable {

    public static final String NETWORK_PARAM_FILE = "network/parameters.txt";
    public static final String PC_LIST_FILE = "network/pcNames.txt";

    public FlowProProperties props;

    public final double h = 1e-8;
//    public static final double R_eps = 1e-2;  // tolerance pro hodnotu hustoty

    public final boolean continueComputation;
    public final double residuum;
    public final double tEnd;
    public final int steps;         // pocet casovych kroku
    public final int saveRate;
    public final boolean animation;
    public final boolean curvedBoundary;

    public double cfl;        // max CFL cislo
    public boolean varyCFL;
    public int order;         // rad metody v prostoru
    public final int timeOrder;         // rad metody v case
    public final int nThreads;      // pocet vlaken
    public final int newtonIters;   // pocet vnitrnich iteraci    
    public final double newtonIterTol;
    public final double meshScale; // mesh scale

    // transformation object
    public DomainTransformationObject domainTransformationObject;

    // integration rules
    public int volumeQuardatureOrder;
    public int faceQuardatureOrder;

    // solver type
    public String linearSolver;
    public String preconditioner;
    public double iterativeSolverTol;

    // single time / dual time
    public String spatialMethod;
    public String timeMethod;
    public String parallelSolverType;
    public String parallelPreconditioner;

    // dynamics parameters
    public boolean movingMesh;
    public String movementType;
    public double bodyLength; // length of body in z direction

    // dynamics model
    public String dynamicsModel;
    public String meshDeformationType;

    /* parallel mode */
    public final boolean parallelMode;
    public final int overlap;
    public final int schwarzIters;
    public final double schwarzTol;
    public final int slavePort;
    public final int fetcherPort;
    public final String pcFilterFile;
    public final String publicKeyFile;

    // solution monitor
    public String solutionMonitor;
    public boolean solutionMonitorOn;
    public boolean solutionAverage; // solution average in elementData

    // external field
    public boolean externalField;

    public Parameters(String parameterFilePath, boolean parallelMode, URL[] jarURLList) throws IOException {
        try {
            this.parallelMode = parallelMode;
            props = new FlowProProperties();
            props.load(new FileInputStream(parameterFilePath));

            continueComputation = props.getBoolean("continueComputation");
            steps = props.getInt("steps");

            if (props.containsKey("residuum")) {
                residuum = props.getDouble("residuum");
            } else {
                residuum = -1.0;
            }

            if (props.containsKey("meshScale")) {
                meshScale = props.getDouble("meshScale");
            } else {
                meshScale = 1;
            }

            if (props.containsKey("tEnd")) {
                tEnd = props.getDouble("tEnd");
            } else {
                tEnd = Double.MAX_VALUE;
            }

            if (props.containsKey("saveRate")) {
                saveRate = props.getInt("saveRate");
            } else {
                saveRate = Integer.MAX_VALUE;
            }

            if (props.containsKey("animation")) {
                animation = props.getBoolean("animation");
            } else {
                animation = false;
            }

            if (props.containsKey("spatialMethod")) {
                spatialMethod = props.getString("spatialMethod");
            } else {
                flowpro.core.Mesh.SpatialMethodType.help();
                throw new IOException("parameter spatialMethod must by defined");
            }

            if (props.containsKey("order")) {
                order = props.getInt("order");
            } else {
                order = 1;
            }

            if (props.containsKey("timeOrder")) {
                timeOrder = props.getInt("timeOrder");
            } else {
                timeOrder = 1;
            }

            if (props.containsKey("timeMethod")) {
                timeMethod = props.getString("timeMethod");
            } else {
                flowpro.core.element.Element.TimeIntegrationElementType.help();
                throw new IOException("parameter timeMethod must by defined");
            }

            volumeQuardatureOrder = order;
            if (props.containsKey("volumeQuardatureOrder")) {
                volumeQuardatureOrder = props.getInt("volumeQuardatureOrder");
            }

            faceQuardatureOrder = order;
            if (props.containsKey("faceQuardatureOrder")) {
                faceQuardatureOrder = props.getInt("faceQuardatureOrder");
            }

            cfl = props.getDouble("CFL");
            if (cfl == -1) {
                varyCFL = true;
                cfl = 1;
            }

            nThreads = props.getInt("threads");
            newtonIters = props.getInt("newtonIters");

            if (props.containsKey("newtonIterTol")) {
                newtonIterTol = props.getDouble("newtonIterTol");
            } else {
                newtonIterTol = 1e-4;
            }

            // domain transformation object
            domainTransformationObject = null;
            if (props.containsKey("domainTransformationObject")) {
                domainTransformationObject = (new DomainTransformationObjectFactory()).getDomainTransformationObject(parameterFilePath, jarURLList);
            }

            // dynamics parameters
            movingMesh = false;
            if (props.containsKey("movingMesh")) {
                movingMesh = props.getBoolean("movingMesh");
            }

            // iterative solver setting
            if (props.containsKey("linearSolver")) {
                linearSolver = props.getString("linearSolver");
            } else {
                flowpro.core.LinearSolvers.LinearSolver.LinearSolverType.help();
                throw new IOException("parameter linearSolver must by defined");
            }

            if (props.containsKey("preconditioner")) {
                preconditioner = props.getString("preconditioner");
            } else {
                flowpro.core.LinearSolvers.preconditioners.Preconditioner.PreconditionerType.help();
                throw new IOException("parameter preconditioner must by defined");
            }

            if (props.containsKey("iterativeSolverTol")) {
                iterativeSolverTol = props.getDouble("iterativeSolverTol");
            } else {
                iterativeSolverTol = 1e-2;
            }

            if (parallelMode) {
                if (props.containsKey("parallelSolverType")) {
                    parallelSolverType = props.getString("parallelSolverType");
                } else {
                    flowpro.core.solver.MasterSolver.MasterSolverType.help();
                    throw new IOException("parameter parallelSolverType must by defined");
                }

                if (props.containsKey("parallelPreconditioner")) {
                    parallelPreconditioner = props.getString("parallelPreconditioner");
                } else {
                    flowpro.core.DistributedLinearSolver.preconditioner.ParallelPreconditioner.ParallelPreconditionerType.help();
                    throw new IOException("parameter parallelPreconditioner must by defined");
                }
            }

            if (props.containsKey("curvedBoundary")) {
                curvedBoundary = props.getBoolean("curvedBoundary");
            } else {
                curvedBoundary = false;
            }

            if (props.containsKey("dynamicsModel")) {
                dynamicsModel = props.getString("dynamicsModel");
            } else {
                dynamicsModel = "none";
            }

            if (props.containsKey("meshDeformationType")) {
                meshDeformationType = props.getString("meshDeformationType");
            } else {
                meshDeformationType = "none";
            }

            /* parallel mode */
            if (parallelMode) {
                overlap = props.getInt("overlap");

                if (props.containsKey("schwarzIters")) {
                    schwarzIters = props.getInt("schwarzIters");
                } else {
                    schwarzIters = 1;
                }

                if (props.containsKey("schwarzTol")) {
                    schwarzTol = props.getDouble("schwarzTol");
                } else {
                    schwarzTol = 1e-6;
                }
            } else {
                overlap = -1;
                schwarzIters = -1;
                schwarzTol = 0.0;
            }

            solutionMonitorOn = false;
            if (props.containsKey("solutionMonitor")) {
                solutionMonitor = props.getString("solutionMonitor");
                solutionMonitorOn = true;
            } else {
                solutionMonitor = null;
            }

            solutionAverage = false;
            if (props.containsKey("solutionAverage")) {
                solutionAverage = props.getBoolean("solutionAverage");
            }

        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage());
        }

        try {
            if (parallelMode) {
                FlowProProperties netProps = new FlowProProperties();
                netProps.load(new FileInputStream(NETWORK_PARAM_FILE));
                slavePort = netProps.getInt("port");
                fetcherPort = netProps.getInt("fetcherPort");
                pcFilterFile = netProps.getString("pcFilterFile");
                publicKeyFile = netProps.getString("publicKeyFile");
            } else {
                slavePort = -1;
                fetcherPort = -1;
                pcFilterFile = null;
                publicKeyFile = null;
            }
        } catch (IOException ex) {
            throw new IOException("file " + NETWORK_PARAM_FILE
                    + " has a wrong format: " + ex.getMessage());
        }
    }
}
