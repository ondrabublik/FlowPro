package flowpro.core;

import flowpro.api.Dynamics;
import flowpro.api.Mat;
import flowpro.api.Equation;
import flowpro.api.SolutionMonitor;
import flowpro.core.curvedBoundary.CurvedBoundary;
import flowpro.core.quadrature.QuadratureCentral;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.parallel.Domain;
import flowpro.core.elementType.ElementType;
import flowpro.core.meshDeformation.*;
import litempi.MPIException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import flowpro.core.solver.MasterSolver;
import flowpro.core.solver.SlaveSolver;
import static flowpro.core.elementType.ElementType.firstDigit;
import static flowpro.core.elementType.ElementType.firstDigit;
import static flowpro.core.elementType.ElementType.firstDigit;
import static flowpro.core.elementType.ElementType.firstDigit;

/**
 *
 * @author ales
 */
public class FlowProMain {

    private static URL[] jarURLList;
    private static final Logger LOG = LoggerFactory.getLogger(FlowProMain.class);

    public static final String ARG_FILE_NAME = "args.txt";
    public static final String PARAMETER_FILE_NAME = "parameters.txt";
    public static final String STATE_FILE_NAME = "state.txt";
    public static final String REF_VALUE_FILE_NAME = "referenceValues.txt";

    public final String meshPath;
    public final String simulationPath;
    private final Object lock;

    public FlowProMain() throws IOException {
        lock = new Object();

        File f = new File(ARG_FILE_NAME);
        if (!f.exists()) {
            LOG.warn("File args.txt not exists. Default args.txt file was create.");
            f.createNewFile();
            PrintWriter writer = new PrintWriter(ARG_FILE_NAME, "UTF-8");
            writer.print("examples/GAMM default");
            writer.close();
        }

        BufferedReader reader;
        reader = new BufferedReader(new FileReader(ARG_FILE_NAME));

        String line;
        String[] args;
        String geometryName, simulationName;
        if ((line = reader.readLine()) != null) {
            args = line.split(" ");
            if (args.length != 2) {
                throw new IOException("file " + ARG_FILE_NAME
                        + " must contain one line with exactly two arguments");
            }
            geometryName = args[0];
            simulationName = args[1];

            LOG.info("geometry: {}, simulation: {}", geometryName, simulationName);
        } else {
            throw new IOException("file " + ARG_FILE_NAME + " is empty");
        }

        meshPath = "simulations/" + geometryName + "/mesh/";
        simulationPath = "simulations/" + geometryName + "/" + simulationName + "/";
    }

    public Object getLock() {
        return lock;
    }

    public static void main(String args[]) throws InterruptedException, IOException {
        System.out.println();
        System.out.println();
        System.out.println("FFFFF  L       OOOO   W    W    W  PPPPP   RRRRR    OOOO ");
        System.out.println("F      L      O    O  W   W W   W  P    P  R    R  O    O");
        System.out.println("FFFFF  L      O    O   W  W W  W   PPPPP   RRRRR   O    O");
        System.out.println("F      L      O    O   W W   W W   P       R   R   O    O");
        System.out.println("F      LLLLL   OOOO     W     W    P       R    R   OOOO ");
        System.out.println();
        System.out.println("---------------------------------------------------------");
        System.out.println("Welcome to the CFD opensource software FlowPro!          ");
        System.out.println("---------------------------------------------------------");
        System.out.println();

        LOG.info("starting FlowPro...");
        jarURLList = getJarURLList("modules");
        try {
            if (args == null || args.length == 0) {
                throw new IllegalArgumentException("no argument was specified");
            }

            FlowProMain dgfem;
            MasterSolver solver;
            SlaveSolver slave;
            Solution solution;
            switch (args[0].toLowerCase()) {
                case "master":
                    dgfem = new FlowProMain();
                    solver = dgfem.solverFactory(false, Integer.valueOf(args[1]));
                    solution = solver.solve();
                    solver.saveData(solution);
                    break;

                case "remote":
                    dgfem = new FlowProMain();
                    Object lock = dgfem.getLock();
                    Proxy proxy = new Proxy(args[1], 6666, args[2], dgfem.simulationPath, lock);
                    Thread proxyThread = new Thread(proxy);
                    proxyThread.start();
                    solver = dgfem.solverFactory(false, 0);
                    solution = solver.solve();
                    solver.saveData(solution);
                    synchronized (lock) {
                        proxy.stop();
                        lock.notify();
                    }
                    break;

                case "slave":
                    if (args.length < 4) {
                        throw new IllegalArgumentException("missing arguments after slave");
                    }
                    String masterIP = args[1];
                    int masterPort = Integer.parseInt(args[2]);
                    String parallelSolverType = args[3];
                    try {
                        slave = SlaveSolver.factory(parallelSolverType, masterIP, masterPort);
                        LOG.info("mesh was received and initialised");
                        slave.solve();
                    } catch (IOException | MPIException ex) {
                        LOG.error("", ex);
                    }
                    break;
//                    proxyThread.join();

                case "postprocessing": // plot data
                    dgfem = new FlowProMain();
                    ResultsPlot plot = new ResultsPlot(dgfem.meshPath, dgfem.simulationPath, args, jarURLList);
                    plot.generateResults();
                    LOG.info("results were generated into " + dgfem.simulationPath + "output directory");
                    break;

                case "optimalisationexport":
                    dgfem = new FlowProMain();
                    solver = dgfem.solverFactory(true, 0);
                    new OptimisationToolExport(solver, dgfem.simulationPath, args[1].toLowerCase(), jarURLList).export();
                    LOG.info("Optimalisation arrays were exported..");
                    break;

                case "testdynamicmodel":
                    dgfem = new FlowProMain();
                    solver = dgfem.solverFactory(false, 0);
                    solver.testDynamic(Double.valueOf(args[1]));
                    break;

                case "getparameters": // get all solver parameters
                    try {
                        dgfem = new FlowProMain();
                        System.out.println();
                        System.out.println("FlowPro parameters:");
                        Class par = new Parameters(dgfem.simulationPath + PARAMETER_FILE_NAME, false, jarURLList).getClass();
                        Field[] fields = par.getFields();
                        for (Field field : fields) {
                            System.out.println(field.getName() + ": " + field.getGenericType());
                        }
                    } catch (Exception ex) {
                        LOG.error(ex.getMessage());
                    }
                    break;

                default:
                    throw new IllegalArgumentException("unknown argument " + args[0]);
            }
        } catch (IllegalArgumentException ex) {
            LOG.error("input arguments have a wrong format: {}", ex.getMessage(), ex);
            System.exit(1);
        } catch (MPIException | IOException ex) {
            LOG.error("{}", ex.getMessage());
        }
    }

    private MasterSolver solverFactory(boolean optimalisation, int nDomains) throws IOException {
        boolean parallelMode = false;
        if (nDomains > 0) {
            parallelMode = true;
        }
        Equation eqn = (new EquationFactory()).getEquation(simulationPath + PARAMETER_FILE_NAME, jarURLList);   // read physical parameters
        Parameters par = new Parameters(simulationPath + PARAMETER_FILE_NAME, parallelMode, jarURLList); // read numerical parameters                    

        LOG.info("loading data...");
        // load matrices defining the mesh
        double[][] PXY = Mat.loadDoubleMatrix(meshPath + "vertices.txt"); // mesh vertices coordinates
        if (par.meshScale != 1) {
            for (int i = 0; i < PXY.length; i++) {
                for (int j = 0; j < PXY[i].length; j++) {
                    PXY[i][j] *= par.meshScale;
                }
            }
            LOG.info("Mesh was scaled with parameter " + par.meshScale);
        }

        double[][] UXY = null;
        if (par.movingMesh && (par.continueComputation || optimalisation)) {
            try {
                UXY = Mat.loadDoubleMatrix(simulationPath + "UXY.txt"); // mesh vertices velocity
                LOG.info("Mesh velocities were loades!");
            } catch (FileNotFoundException ex) {
                LOG.warn("Mesh velocities file UXY.txt was not found. Velocities were set to zero!");
            }
        }

        int[] elemsType = Mat.loadIntArray(meshPath + "elementType.txt");   // element type
        int[][] TP = Mat.loadIntMatrix(meshPath + "elements.txt"); // inexes of points defining element
        LOG.info("Mesh info: " + PXY.length + "-vertices, " + TP.length + "-elements.");

        int[][] TT;
        try {
            TT = Mat.loadIntMatrix(meshPath + "neighbors.txt");  // indexes of element neigborhoods
        } catch (FileNotFoundException ex) {
            LOG.info("creating neighbor matrix");
            int[][] boundary = Mat.loadIntMatrix(meshPath + "boundaryType.txt");
            TT = createNeighborMatrix(elemsType, TP, boundary);
            Mat.save(TT, meshPath + "neighbors.txt");
        }
        int nElems = TP.length;  // number of elements in the mesh

        int[] elemsOrder;
        try {
            elemsOrder = Mat.loadIntArray(simulationPath + "order.txt");
            LOG.info("reading local order of spatial accuracy from file " + "order.txt");
        } catch (FileNotFoundException ex) {
            if (par.order < 1) {
                throw new IOException("neither global nor local order of spatial accuracy defined, "
                        + " either define variable order in file " + simulationPath + PARAMETER_FILE_NAME
                        + " or create file " + "order.txt" + " in simulation path");
            }
            elemsOrder = new int[nElems];
            LOG.warn("file " + simulationPath + "order.txt" + " not found"
                    + ", setting global order of spatial accuracy to " + par.order);
            Arrays.fill(elemsOrder, par.order);
        }

        // curved boundary computation
        FaceCurvature[] fCurv;
        if (par.curvedBoundary) {
            fCurv = CurvedBoundary.modifyMesh(elemsType, PXY, TP, TT);
            CurvedBoundary.saveMesh(meshPath, elemsType, TP, TT);
        } else {
            fCurv = new FaceCurvature[nElems];
            //elemsType = firstDigit(elemsType);
        }

        // domain partition for parallel computing
        int[] elem2DomMap;
        if (parallelMode) {
            try {
                elem2DomMap = Mat.loadIntArray(simulationPath + "part.txt");
            } catch (FileNotFoundException ex) {
                LOG.info("generating domain partition");
                elem2DomMap = generateMeshPartition(elemsType, PXY, TP, nDomains);
                Mat.save(elem2DomMap, simulationPath + "part.txt");
            }
        } else {
            elem2DomMap = new int[nElems];
        }

        // wall distance for turbulence models
        double[] wallDistance;
        try {
            wallDistance = Mat.loadDoubleArray(meshPath + "wallDistance.txt");
        } catch (FileNotFoundException ex) {
            LOG.info("generating wall distance function");
            wallDistance = generateWallDistanceFunction(eqn, elemsType, PXY, TP, TT);
            Mat.save(wallDistance, meshPath + "wallDistance.txt");
        }

        // body boundaries for ALE computation
        int[][] TEale;
        if (par.movingMesh) {
            LOG.info("moving mesh computation - ALE description");
            try {
                TEale = Mat.loadIntMatrix(meshPath + "neighborsALE.txt");
            } catch (FileNotFoundException ex) {
                LOG.info("creating neighborALE matrix");
                int[][] boundaryALE = Mat.loadIntMatrix(meshPath + "boundaryTypeALE.txt");
                TEale = createNeighborALEMatrix(elemsType, TP, TT, boundaryALE);
                Mat.save(TEale, meshPath + "neighborsALE.txt");
            }
        } else {
            LOG.info("static mesh computation - Eulerian description");
            TEale = Mat.allocSameIntMatrix(TT);
        }

        // shifting for periodicity boundary condition
        int[][] TEshift;
        double[][] shift = null;
        try {
            TEshift = Mat.loadIntMatrix(meshPath + "TEshift.txt");
            shift = Mat.loadDoubleMatrix(meshPath + "shift.txt");
            if (par.meshScale != 1) {
                for (int i = 0; i < shift.length; i++) {
                    for (int j = 0; j < shift[i].length; j++) {
                        shift[i][j] *= par.meshScale;
                    }
                }
            }
            LOG.info("periodic boundary found");
        } catch (FileNotFoundException ex) {
            TEshift = Mat.allocSameIntMatrix(TT);
            LOG.info("no periodic boundary found");
        }

        int nDoms = Mat.max(elem2DomMap) + 1;

        printInfo(eqn, par);

        // loading external field
        par.externalField = false;
        double[][] externalField = null;
        try {
            externalField = Mat.loadDoubleMatrix(meshPath + "externalField.txt");
            par.externalField = true;
            LOG.info(" external field found");
        } catch (FileNotFoundException ex) {
            LOG.info(" external field not found ");
        }

        // load initial condition
        double initW[][];
        State state = new State(simulationPath + STATE_FILE_NAME, par.order);
        if (optimalisation) {
            initW = Mat.loadDoubleMatrix(simulationPath + "We.txt");
        } else if (par.continueComputation) {
            state.load();
            if (state.hasOrderChanged()) {
                initW = Mat.loadDoubleMatrix(simulationPath + "W.txt");
            } else {
                try {
                    initW = Mat.loadDoubleMatrix(simulationPath + "We.txt");
                } catch (FileNotFoundException ex) {
                    LOG.warn("file " + simulationPath + "We.txt" + " not found, "
                            + "reading " + simulationPath + "W.txt" + " instead");
                    initW = Mat.loadDoubleMatrix(simulationPath + "W.txt");
                }
            }
        } else {
            try {
                initW = Mat.loadDoubleMatrix(simulationPath + "initW.txt");
                LOG.info("init conditions were loaded from file initW.txt ");
            } catch (IOException e) {
                try {
                    System.out.println("Trying run script " + simulationPath + "initScript.js");
                    initW = (new InitConditionScript(PXY, TP, eqn.nEqs(), simulationPath + "initScript.js")).getInitW();
                } catch (Exception ex) {
                    System.out.println("Script initScript.js cannot be expressed because: " + ex);
                    double[] w = eqn.constInitCondition();
                    initW = replicate(w, nElems, eqn.nEqs());
                }
            }
        }

        // create solution monitor
        SolutionMonitor solMonitor = null;
        if (par.solutionMonitorOn) {
            solMonitor = (new SolutionMonitorFactory()).getSolutionMonitor(simulationPath + PARAMETER_FILE_NAME, jarURLList);
        }

        // create file for residuum
        if (!par.continueComputation) {
            try {
                PrintWriter writer = new PrintWriter(simulationPath + "residuum.txt");
                writer.print("");
                writer.close();
            } catch (IOException ioe) {
                LOG.error("error while creating a new empty file residuum.txt");
            }
        }

        // load Gaussian quadrature
        QuadratureCentral qRules = new QuadratureCentral();

        LOG.info("initialising mesh...");
        // initialising ALE objects
        Deformation dfm = new DeformationFactory().getDeformation(par, eqn, TEale);
        Dynamics dyn = null;
        if (par.movingMesh) {
            dyn = (new DynamicsFactory()).getDynamicsModel(simulationPath + PARAMETER_FILE_NAME, jarURLList, dfm.nBodies, simulationPath, meshPath, eqn);
            dfm.setCenters(dyn.getCenter());
            dfm.calculateBlendingFunctions(PXY, TP, TT, TEale, elemsType, meshPath);
        }

        Domain domain = new Domain(elemsOrder, elemsType, TT, TP, TEale, TEshift, fCurv, initW, elem2DomMap, nDoms, par.overlap, PXY.length);
        Mesh[] mesh = new Mesh[nDoms];
        for (int d = 0; d < nDoms; ++d) {
            LOG.info("domain " + d + " from " + (nDoms - 1) + ":");
            Domain.Subdomain subdom = domain.getSubdomain(d);
            mesh[d] = new Mesh(eqn, dfm, par, solMonitor, qRules, PXY, UXY, subdom.elemsOrder, wallDistance, externalField, subdom.elemsType,
                    subdom.TP, subdom.TT, subdom.TEale, subdom.TEshift, shift, subdom.fCurv, subdom.initW, subdom);
            if (!parallelMode) {
                mesh[d].init();
            }
        }

        return MasterSolver.factory(simulationPath, mesh, dyn, eqn, par, state, domain, lock);
    }

    public static String millisecsToTime(long nanoseconds) {
        long seconds = nanoseconds / 1000;

        int hours = (int) seconds / 3600;
        int minutes = (int) seconds / 60 - hours * 60;
        seconds = (int) seconds - minutes * 60 - hours * 3600;

        if (hours == 0 && minutes == 0) {
            return String.format("%02ds", seconds);
        } else if (hours == 0) {
            return String.format("%02dm %02ds", minutes, seconds);
        } else {
            return String.format("%02dh %02dm %02ds", hours, minutes, seconds);
        }
    }

    private static void printInfo(Equation eqn, Parameters par) {
        LOG.info("CFL: {}, order: {}, threads: {}", par.cfl, par.order, par.nThreads);
        StringBuilder strBuilder = new StringBuilder();
        strBuilder.append("steps: ").append(par.steps);
        if (par.residuum > 0) {
            strBuilder.append(", ").append("residuum: ").append(par.residuum);
        }
        if (par.tEnd < 1e8) {
            strBuilder.append(", ").append("tEnd: ").append(par.tEnd);
        }
        LOG.info(strBuilder.toString());
    }

    private static double[][] replicate(double[] initW, int nElems, int nEqs) {
        double[][] W = new double[nElems][];
        for (int i = 0; i < nElems; ++i) {
            W[i] = new double[nEqs];
            for (int j = 0; j < nEqs; ++j) {
                W[i][j] = initW[j];
            }
        }
        return W;
    }

    public int[][] createNeighborMatrix(int[] type, int[][] TP, int[][] boundary) {

        // neighbors finder
        HashMap<String, Object> map = new HashMap<>();

        int nRows = TP.length;
        // creating elementTypes
        ElementType[] elemType = new ElementType[nRows];
        for (int i = 0; i < nRows; i++) {
            elemType[i] = ElementType.elementTypeFactory(type[i], 0, 0, 0);
        }

        // saving face indexes
        for (int i = 0; i < nRows; i++) {
            int nFaces = elemType[i].numberOfEdges();
            for (int k = 0; k < nFaces; k++) {
                int[] faceIndexes = elemType[i].getFaceIndexes(k);
                for (int j = 0; j < faceIndexes.length; j++) {
                    faceIndexes[j] = TP[i][faceIndexes[j]];
                }
                Arrays.sort(faceIndexes);
                String key = "";
                for (int j = 0; j < faceIndexes.length; j++) {
                    key = key + Integer.toString(faceIndexes[j]) + "$";
                }
                FaceIndexes face = (FaceIndexes) map.get(key);
                if (face == null) {
                    face = new FaceIndexes(i);
                    map.put(key, face);
                } else {
                    face.i2 = i;
                }
            }
        }

        // loading face indexes
        int[][] TT = new int[nRows][];
        for (int i = 0; i < nRows; i++) {
            int nFaces = elemType[i].numberOfEdges();
            TT[i] = new int[nFaces];
            Arrays.fill(TT[i], -1000);
            for (int k = 0; k < nFaces; k++) {
                int[] faceIndexes = elemType[i].getFaceIndexes(k);
                for (int j = 0; j < faceIndexes.length; j++) {
                    faceIndexes[j] = TP[i][faceIndexes[j]];
                }
                Arrays.sort(faceIndexes);
                String key = "";
                for (int j = 0; j < faceIndexes.length; j++) {
                    key = key + Integer.toString(faceIndexes[j]) + "$";
                }
                FaceIndexes face = (FaceIndexes) map.get(key);
                if (face.i1 == i) {
                    if (face.i2 != -1000) {
                        TT[i][k] = face.i2;
                    }
                } else {
                    TT[i][k] = face.i1;
                }
            }
        }

        // boundary
        HashMap<String, Integer> boundMap = new HashMap<>();
        int nBoundRows = boundary.length;
        for (int i = 0; i < nBoundRows; i++) {
            int[] faceIndexes = new int[boundary[i].length - 1];
            for (int j = 0; j < boundary[i].length - 1; j++) {
                faceIndexes[j] = boundary[i][j + 1];
            }
            Arrays.sort(faceIndexes);
            String key = "";
            for (int j = 0; j < faceIndexes.length; j++) {
                key = key + Integer.toString(faceIndexes[j]) + "$";
            }
            boundMap.put(key, boundary[i][0]);
        }

        for (int i = 0; i < nRows; i++) {
            for (int k = 0; k < TT[i].length; k++) {
                if (TT[i][k] == -1000) {
                    int[] faceIndexes = elemType[i].getFaceIndexes(k);
                    for (int j = 0; j < faceIndexes.length; j++) {
                        faceIndexes[j] = TP[i][faceIndexes[j]];
                    }
                    Arrays.sort(faceIndexes);
                    String key = "";
                    for (int j = 0; j < faceIndexes.length; j++) {
                        key = key + Integer.toString(faceIndexes[j]) + "$";
                    }
                    Integer index = boundMap.get(key);
                    if (index != null) {
                        TT[i][k] = index;
                    } else {
                        TT[i][k] = -1;
                    }
                }
            }
        }

        return TT;
    }

    public int[][] createNeighborALEMatrix(int[] type, int[][] TP, int[][] TT, int[][] boundary) {

        int nRows = TP.length;
        // creating elementTypes
        ElementType[] elemType = new ElementType[nRows];
        for (int i = 0; i < nRows; i++) {
            elemType[i] = ElementType.elementTypeFactory(type[i], 0, 0, 0);
        }

        int[][] TEale = new int[TT.length][];

        HashMap<String, Integer> boundMap = new HashMap<>();
        int nBoundRows = boundary.length;
        for (int i = 0; i < nBoundRows; i++) {
            int[] faceIndexes = new int[boundary[i].length - 1];
            for (int j = 0; j < boundary[i].length - 1; j++) {
                faceIndexes[j] = boundary[i][j + 1];
            }
            Arrays.sort(faceIndexes);
            String key = "";
            for (int j = 0; j < faceIndexes.length; j++) {
                key = key + Integer.toString(faceIndexes[j]) + "$";
            }
            boundMap.put(key, boundary[i][0]);
        }

        for (int i = 0; i < TT.length; i++) {
            TEale[i] = new int[TT[i].length];
            for (int k = 0; k < TT[i].length; k++) {
                //if (TT[i][k] < 0) {
                int[] faceIndexes = elemType[i].getFaceIndexes(k);
                for (int j = 0; j < faceIndexes.length; j++) {
                    faceIndexes[j] = TP[i][faceIndexes[j]];
                }
                Arrays.sort(faceIndexes);
                String key = "";
                for (int j = 0; j < faceIndexes.length; j++) {
                    key = key + Integer.toString(faceIndexes[j]) + "$";
                }
                Integer index = boundMap.get(key);
                if (index != null) {
                    TEale[i][k] = index;
                } else {
                    TEale[i][k] = 0;
                }
                //}
            }
        }
        return TEale;
    }

    static class FaceIndexes {

        int i1, i2;

        FaceIndexes(int i) {
            i1 = i;
            i2 = -1000;
        }
    }

    public double[] generateWallDistanceFunction(Equation eqn, int[] type, double[][] PXY, int[][] TP, int[][] TT) {
        int nRows = TP.length;
        // creating elementTypes
        ElementType[] elemType = new ElementType[nRows];
        for (int i = 0; i < nRows; i++) {
            elemType[i] = ElementType.elementTypeFactory(type[i], 0, 0, 0);
        }

        int nPoints = PXY.length;
        int nElem = TT.length;
        double[] wallDistance = new double[nPoints];
        Arrays.fill(wallDistance, 1e30);
        for (int i = 0; i < nElem; i++) {
            for (int k = 0; k < TT[i].length; k++) {
                if (eqn.isIPFace(TT[i][k])) {
                    int[] faceIndexes = elemType[i].getFaceIndexes(k);
                    for (int j = 0; j < faceIndexes.length; j++) {
                        wallDistance[TP[i][faceIndexes[j]]] = 0;
                    }
                }
            }
        }

        int iter = 1;
        double sum = 0;
        while (true) {
            for (int i = 0; i < nElem; i++) {
                for (int j = 0; j < TP[i].length; j++) {
                    if (wallDistance[TP[i][j]] > 0) {
                        for (int k = 0; k < TP[i].length; k++) {
                            if (j != k) {
                                double L = 0;
                                for (int d = 0; d < PXY[0].length; d++) {
                                    L += (PXY[TP[i][j]][d] - PXY[TP[i][k]][d]) * (PXY[TP[i][j]][d] - PXY[TP[i][k]][d]);
                                }
                                L = Math.sqrt(L);
                                double y = wallDistance[TP[i][k]] + L;
                                if (y < wallDistance[TP[i][j]]) {
                                    wallDistance[TP[i][j]] = y;
                                }
                            }
                        }
                    }
                }
            }
            double sumOld = sum;
            sum = 0;
            for (int i = 0; i < nPoints; i++) {
                sum += wallDistance[i];
            }

            if (sum < 1e30 && Math.abs(sumOld - sum) == 0) {
                System.out.println("wall distance generated at " + iter + " th iteration");
                break;
            }

            iter++;
            if (iter > 1e4) {
                System.out.println("wall distance error, probably wall boundary condition not defined");
                break;
            }
        }
        return wallDistance;
    }

    public int[] generateMeshPartition(int[] elemsType, double[][] PXY, int[][] TP, int nDoms) {
        int[] part = new int[TP.length];
        int dim = PXY[0].length;

        // finding largest dimension of mesh
        double[] dist = new double[dim];
        int[] index = new int[dim];
        for (int d = 0; d < dim; d++) {
            double min = 1e10;
            double max = -1e10;
            for (int i = 0; i < PXY.length; i++) {
                if (min < PXY[i][d]) {
                    min = PXY[i][d];
                }
                if (max > PXY[i][d]) {
                    max = PXY[i][d];
                }
            }
            dist[d] = max - min;
            index[d] = d;
        }
        Mat.quicksort(dist, index, 0, dim - 1);
        int largestDim = index[dim - 1];

        // generating elements centers and sorting
        double[] elemCenter = new double[TP.length];
        int[] indexes = new int[TP.length];
        for (int i = 0; i < TP.length; i++) {
            for (int j = 0; j < TP[i].length; j++) {
                elemCenter[i] += PXY[TP[i][j]][largestDim];
            }
            elemCenter[i] /= TP[i].length;
            indexes[i] = i;
        }
        Mat.quicksort(elemCenter, indexes, 0, elemCenter.length - 1);

        // make partitioning
        int mapInd = 0;
        int s = 0;
        for (int i = 0; i < TP.length; i++) {
            part[indexes[i]] = mapInd;
            s = s + 1;
            if (s > TP.length / nDoms) {
                mapInd = mapInd + 1;
                s = 0;
            }
        }
        return part;
    }

    public static URL[] getJarURLList(String s) throws IOException {
        URL[] u = null;
        try {
            File currentDir = new File(s); // current directory
            ArrayList<URL> URLs = new ArrayList();
            System.out.print("Found libraries: ");
            addDirectoryContents(currentDir, URLs);
            System.out.println();
            u = new URL[URLs.size()];
            URLs.toArray(u);
        } catch (Exception e) {
            System.out.println(e);
        }
        return u;
    }

    public static void addDirectoryContents(File dir, ArrayList<URL> URLs) {
        try {
            File[] files = dir.listFiles();
            for (File file : files) {
                if (file.isDirectory()) {
                    addDirectoryContents(file, URLs);
                } else {
                    URLs.add(file.toURI().toURL());
                    System.out.print(file.getName() + ", ");
                }
            }
        } catch (Exception e) {
            System.out.println("Error " + e);
        }
    }
}
