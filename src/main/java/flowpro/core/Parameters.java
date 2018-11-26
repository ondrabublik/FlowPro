package flowpro.core;

import flowpro.api.FlowProProperties;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.Serializable;

/**
 *
 * @author ales
 */
public class Parameters implements Serializable {

    public static final String NETWORK_PARAM_FILE = "network/parameters.txt";

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
    public double cflLTS;
    public boolean varyCFL;
    public final int order;         // rad metody v prostoru
    public int orderInTime; // rad metody v case
    public final int nThreads;      // pocet vlaken
    public final int newtonIters;   // pocet vnitrnich iteraci    
    public final double newtonIterTol;
    public final double penalty;    // interior penalty constant
    public final double beta0;  // direct discontinuous constant 
    public final double meshScale; // mesh scale
    public boolean useJacobiMatrix;
    
    // solver type
    public String linearSolver;
    public String preconditioner;
    public double iterativeSolverTol;

    // single time / dual time
    public String localSolverType;
    public String parallelSolverType;
    public boolean isExplicit;
    public String timeMethod;
    public double[] coeffsPhys;
    public double[] coeffsDual;

    // dynamics parameters
    public boolean movingMesh;
    public String movementType;
    public double bodyLength; // length of body in z direction
    
    // artificial damping
    public final double dampTol;   // tolerance pro pridavnou viskozitu
    public final double dampConst; // konstantni pridavna viskozita 

    // Finite volume method limiter
    public final String FVMlimiter;
    
    // dynamics model
    public String dynamicsModel;
    public String meshDeformationType;
    
    /* parallel mode */
    public final boolean parallelMode;
    public final int overlap;
    public final int schwarzIters;
    public final double schwarzTol;
    public final String masterIP;
    public final int masterPort;
    public final int fetcherPort;

    // solution monitor
    public String solutionMonitor;
    public boolean solutionMonitorOn;
    
    // external field
    public boolean externalField;
    
    public Parameters(String parameterFilePath, boolean parallelMode) throws IOException {
        try {
            this.parallelMode = parallelMode;
            FlowProProperties props = new FlowProProperties();
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

            if (props.containsKey("order")) {
                order = props.getInt("order");
            } else {
                order = Integer.MIN_VALUE;
            }
            
            orderInTime = props.getInt("orderInTime");
            cfl = props.getDouble("CFL");
            if(cfl == -1){
                varyCFL = true;
                cfl = 1;
            }
            
            if (props.containsKey("CFLlts")) {
                cflLTS = props.getDouble("CFLlts");
            } else {
                cflLTS = 0.5;
            }
            
            nThreads = props.getInt("threads");
            newtonIters = props.getInt("newtonIters");
            
            if (props.containsKey("newtonIterTol")) {
                newtonIterTol = props.getDouble("newtonIterTol");
            } else {
                newtonIterTol = 1e-4;
            }
            
            if (props.containsKey("penalty")) {
                penalty = props.getDouble("penalty");
            } else {
                penalty = 0;
            }
            
            if (props.containsKey("beta0")) {
                beta0 = props.getDouble("beta0");
            } else {
                beta0 = 2;
            }

            if (props.containsKey("FVMlimiter")) {
                FVMlimiter = props.getString("FVMlimiter");
            } else {
                FVMlimiter = "none";
            }
            
            if (props.containsKey("useJacobiMatrix")) {
                useJacobiMatrix = props.getBoolean("useJacobiMatrix");
            } else {
                useJacobiMatrix = false;
            }
            
            // dynamics parameters
            movingMesh = false;
            if (props.containsKey("movingMesh")) {
                movingMesh = props.getBoolean("movingMesh");
            }

            // damping
            dampTol = props.getDouble("dampTol");
            dampConst = props.getDouble("dampConst");

            // iterative solver setting
            if (props.containsKey("linearSolver")) {
                linearSolver = props.getString("linearSolver");
            } else {
                linearSolver = "gmres";
            }
            if (props.containsKey("preconditioner")) {
                preconditioner = props.getString("preconditioner");
            } else {
                preconditioner = "blockjacobiinversion";
            }
            if (props.containsKey("iterativeSolverTol")) {
                iterativeSolverTol = props.getDouble("iterativeSolverTol");
            } else {
                iterativeSolverTol = 1e-2;
            }

            localSolverType = "localimplicit";
            isExplicit = false;
            if (props.containsKey("localSolverType")) {
                localSolverType = props.getString("localSolverType");
            }
            
            parallelSolverType = "distKSP";
            if (props.containsKey("parallelSolverType")) {
                parallelSolverType = props.getString("parallelSolverType");
            }
            
            if (props.containsKey("timeMethod")) {                
                timeMethod = props.getString("timeMethod");
                if (timeMethod.equals("dualTime")) {
                    coeffsPhys = props.getDoubleArray("coeffsPhys");
                    coeffsDual = props.getDoubleArray("coeffsDual");
                }
            } else {
                timeMethod = "singleTime";
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

        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage());
        }

        try {
            if (parallelMode) {
                FlowProProperties netProps = new FlowProProperties();
                netProps.load(new FileInputStream(NETWORK_PARAM_FILE));
                masterIP = netProps.getString("ip");
                masterPort = netProps.getInt("port");
                fetcherPort = netProps.getInt("fetcherPort");
            } else {
                masterIP = null;
                masterPort = -1;
                fetcherPort = -1;
            }
        } catch (IOException ex) {
            throw new IOException("file " + NETWORK_PARAM_FILE
                    + " has a wrong format: " + ex.getMessage());
        }
    }
}
