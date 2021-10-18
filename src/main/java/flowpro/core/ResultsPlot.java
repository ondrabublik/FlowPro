package flowpro.core;

import flowpro.api.Equation;
import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import static flowpro.core.FlowProMain.PARAMETER_FILE_NAME;
import flowpro.core.basis.Basis;
import flowpro.core.curvedBoundary.CurvedBoundary;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.elementType.ElementType;
import static flowpro.core.elementType.ElementType.firstDigit;
import flowpro.core.transformation.Transformation;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Stack;

/**
 *
 * @author obublik
 */
public class ResultsPlot {

    String meshPath;
    String simulationPath;
    String outputPath;
    String[] args;
    URL[] jarURLList;

    Equation eqn;

    double[][] Wcoef;
    double[][] PXY;
    int[] elemsType;
    FaceCurvature[] fCurv;
    int[][] TP;
    int[][] TT;
    int dim;
    int[] order;

    FlowProProperties props;

    double lRef;
    double meshScale;

    ResultsPlot(String meshPath, String simulationPath, String[] args, URL[] jarURLList) throws IOException {
        this.meshPath = meshPath;
        this.simulationPath = simulationPath;
        this.args = args;
        this.jarURLList = jarURLList;
        outputPath = simulationPath + "/output/";
    }

    void generateResults() throws IOException {

        Parameters par = new Parameters(simulationPath + PARAMETER_FILE_NAME, false, jarURLList); // read numerical parameters    

        // equations model
        Equation eqn = null;
        try {
            eqn = (new EquationFactory()).getEquation(simulationPath + "parameters.txt", jarURLList);   // read physical parameters
        } catch (Exception e) {
            System.out.println("Error in model." + e);
        }

        lRef = 1;
        try {
            double[] refValues = eqn.getReferenceValues();
            lRef = refValues[0];
            System.out.println("Reference length found: " + lRef + " ");
        } catch (Exception e) {
            System.out.println("Reference length not defined!");
        }

        // loading mesh
        int[][] TP = null;
        try {
            TP = Mat.loadIntMatrix(meshPath + "elements.txt"); // inexes of points defining element
        } catch (FileNotFoundException ex) {
            System.out.println("Elements file not found.");
        }

        // result generating
        String format = "txt";
        int iter = -1;
        int precision = 0;
        int smoothingIter = 0;
        String subName = "";

        int nArgs = args.length;
        Stack<String> names = new Stack();
        for (int i = 1; i < nArgs; i++) {
            switch (args[i].charAt(0)) {
                case '-':
                    switch (args[i].charAt(1)) {
                        case 'f': // format of output results
                            format = args[i].substring(2, args[i].length());
                            break;
                        case 'i': // iteration
                            iter = Integer.valueOf(args[i].substring(2, args[i].length()));
                            break;
                        case 'p': // precision
                            precision = Integer.valueOf(args[i].substring(2, args[i].length()));
                            break;
                        case 's': // results smoothing
                            smoothingIter = Integer.valueOf(args[i].substring(2, args[i].length()));
                            break;
                        case 'n': // sub name
                            subName = args[i].substring(2, args[i].length());
                            break;
                    }
                    break;
                default:
                    names.push(args[i]);
                    break;
            }
        }

        props = new FlowProProperties();
        props.load(new FileInputStream(simulationPath + PARAMETER_FILE_NAME));

        order = null;
        try {
            order = Mat.loadIntArray(simulationPath + "order.txt");
        } catch (FileNotFoundException ex) {
            int orderGlob = props.getInt("order");
            order = new int[TP.length];
            for (int i = 0; i < TP.length; i++) {
                order[i] = orderGlob;
            }
        }

        boolean movingMesh = par.movingMesh;

        // mesh points
        double[][] PXY = null;
        try {
            if (iter != -1 && movingMesh) {
                PXY = Mat.loadDoubleMatrix(simulationPath + "animation/vertices" + String.format("%08d", iter) + ".txt"); // mesh vertices coordinates
            } else {
                PXY = Mat.loadDoubleMatrix(meshPath + "vertices.txt"); // mesh vertices coordinates
            }
        } catch (FileNotFoundException ex) {
            System.out.println("Vertices file not found.");
        }
        dim = PXY[0].length;

        meshScale = 1;
        if (props.containsKey("meshScale")) {
            meshScale = props.getDouble("meshScale");
        }

        lRef = 1;
        if (props.containsKey("lRef")) {
            lRef = props.getDouble("lRef");
        }

        boolean[] scaleDims = new boolean[3];
        Arrays.fill(scaleDims, Boolean.TRUE);
        if (props.containsKey("scaleDims")) {
            scaleDims = props.getBooleanArray("scaleDims");
        } else {
            
        }

        if (meshScale != 1 || lRef != 1) {
            for (int i = 0; i < PXY.length; i++) {
                for (int j = 0; j < PXY[i].length; j++) {
                    if (scaleDims[j]) {
                        PXY[i][j] *= par.meshScale / par.lRef;
                    }
                }
            }
        }

        try {
            order = Mat.loadIntArray(simulationPath + "order.txt");
            System.out.println("reading local order of spatial accuracy from file order.txt");
        } catch (FileNotFoundException ex) {
//            if (par.order < 1) {
//                throw new IOException("neither global nor local order of spatial accuracy defined, "
//                        + " either define variable order in file " + simulationPath + PARAMETER_FILE_NAME
//                        + " or create file " + "order.txt" + " in simulation path");
//            }
            order = new int[TP.length];
            if (props.containsKey("order")) {
                Arrays.fill(order, props.getInt("order"));
            } else {
                throw new IOException("neither global nor local order of spatial accuracy defined, "
                        + " either define variable order in file " + simulationPath + PARAMETER_FILE_NAME
                        + " or create file " + "order.txt" + " in simulation path");
            }

            System.out.println("file " + simulationPath + "order.txt not found"
                    + ", setting global order of spatial accuracy to " + order[0]);
        }

        // loading result
        double[][] W = null;
        try {
            switch (precision) {
                case 0:
                    if (iter != -1) {
                        W = Mat.loadDoubleMatrix(simulationPath + "animation/W" + String.format("%08d", iter) + ".txt");
                    } else {
                        W = Mat.loadDoubleMatrix(simulationPath + "W.txt");
                    }
                    break;
                case 1:
                    if (iter != -1) {
                        W = Mat.loadDoubleMatrix(simulationPath + "animation/We" + String.format("%08d", iter) + ".txt");
                    } else {
                        W = Mat.loadDoubleMatrix(simulationPath + "We.txt");
                    }
                    elemsType = Mat.loadIntArray(meshPath + "elementType.txt");

                    break;
                default:
                    if (iter != -1) {
                        W = Mat.loadDoubleMatrix(simulationPath + "animation/We" + String.format("%08d", iter) + ".txt");
                    } else {
                        W = Mat.loadDoubleMatrix(simulationPath + "We.txt");
                    }
                    //PXY = Mat.loadDoubleMatrix(meshPath + "vertices.txt");
                    elemsType = Mat.loadIntArray(meshPath + "elementType.txt");
                    TP = Mat.loadIntMatrix(meshPath + "elements.txt");
                    TT = Mat.loadIntMatrix(meshPath + "neighbors.txt");
                    Wcoef = Mat.loadDoubleMatrix(simulationPath + "We.txt");

                    if (par.curvedBoundary) {
                        fCurv = CurvedBoundary.modifyMesh(elemsType, PXY, TP, TT);
                    } else {
                        fCurv = new FaceCurvature[TP.length];
                        //elemsType = firstDigit(elemsType);
                    }
                    break;
            }
        } catch (FileNotFoundException ex) {
            System.out.println("Results file not found.");
        }

        if (precision < 2) {
            if (precision == 0) { // only copy center value
                System.out.println("Warning: derivatives dW/dxi are supported only with -p option and for order higher then 1!");
            }
            // result generated and mesh interpolation
            int listSize = names.size();
            for (int s = 0; s < listSize; s++) {
                String variableName = names.pop();
                double[] value = eqn.getResults(W[0], new double[eqn.nEqs() * dim], new double[dim], variableName);
                double[][] result;
                if (value.length == 1) {
                    result = new double[PXY.length][1]; // scalar
                } else {
                    result = new double[PXY.length][3]; // vector
                }
                int[] nNeighElem = new int[PXY.length];
                if (precision == 0) { // only copy center value
                    for (int i = 0; i < TP.length; i++) {
                        double[] XCenter = new double[dim];
                        for (int d = 0; d < dim; d++) {
                            for (int j = 0; j < TP[i].length; j++) {
                                XCenter[d] += PXY[TP[i][j]][d];
                            }
                            XCenter[d] /= TP[i].length;
                        }
                        value = eqn.getResults(W[i], new double[eqn.nEqs() * dim], XCenter, variableName);
                        for (int j = 0; j < TP[i].length; j++) {
                            for (int k = 0; k < value.length; k++) {
                                result[TP[i][j]][k] += value[k];
                            }
                            nNeighElem[TP[i][j]]++;
                        }
                    }
                } else { // copy value at given point
                    for (int i = 0; i < TP.length; i++) {
                        double[][] vertices = new double[TP[i].length][dim];
                        for (int j = 0; j < TP[i].length; j++) {
                            System.arraycopy(PXY[TP[i][j]], 0, vertices[j], 0, dim);
                        }

                        ElementType elemType = ElementType.elementTypeFactory(elemsType[i], order[i], par.volumeQuardatureOrder, par.faceQuardatureOrder);
                        Transformation transform = elemType.getVolumeTransformation(vertices, null, par);

                        Basis basis = elemType.getBasis(transform);
                        int nBasis = basis.nBasis;
                        for (int j = 0; j < TP[i].length; j++) {
                            int nEqs = eqn.nEqs();
                            double[] xiCoord = transform.getXi(PXY[TP[i][j]]);
                            // basis derivative transform
                            double[][] dXiBasis = new double[nBasis][dim];
                            for (int m = 0; m < nBasis; m++) {
                                for (int d = 0; d < dim; d++) {
                                    dXiBasis[m][d] = basis.derBasis(m, xiCoord, d);
                                }
                            }
                            double[][] iT = Mat.invert(transform.jacobiMatrix(xiCoord));
                            double[][] dXBasis = new double[nBasis][dim];
                            for (int m = 0; m < nBasis; m++) {
                                for (int r = 0; r < dim; r++) {
                                    for (int t = 0; t < dim; t++) {
                                        dXBasis[m][r] += iT[t][r] * dXiBasis[m][t];
                                    }
                                }
                            }

                            double[] Wpoint = new double[nEqs];
                            double[] dWpoint = new double[dim * nEqs];
                            for (int k = 0; k < Wpoint.length; k++) {
                                for (int m = 0; m < nBasis; m++) {
                                    Wpoint[k] += W[i][k * nBasis + m] * basis.basisFun(m, xiCoord);
                                    for (int d = 0; d < dim; d++) {
                                        dWpoint[nEqs * d + k] += W[i][k * nBasis + m] * dXBasis[m][d];
                                    }
                                }
                            }
                            value = eqn.getResults(Wpoint, dWpoint, vertices[j], variableName);
                            for (int k = 0; k < value.length; k++) {
                                result[TP[i][j]][k] += value[k];
                            }
                            nNeighElem[TP[i][j]]++;
                        }
                    }
                }

                for (int i = 0; i < PXY.length; i++) {
                    for (int d = 0; d < eqn.dim(); d++) {
                        PXY[i][d] *= meshScale;
                    }
                    for (int k = 0; k < result[i].length; k++) {
                        result[i][k] /= nNeighElem[i];
                    }
                }

                // mesh saving
                File directory = new File(outputPath);
                if (!directory.exists()) {
                    directory.mkdir();
                }

                if (smoothingIter > 0) {
                    smoothResults(result, TP, smoothingIter);
                }

                switch (format) {
                    case "txt":
                        Mat.save(result, outputPath + variableName + subName + ".txt");
                        break;
                    case "vtk":
                        int[] elementType = Mat.loadIntArray(meshPath + "elementType" + ".txt");
                        if (s == 0) {
                            save2VTK(result, PXY, TP, elementType, variableName, outputPath + "results" + subName + ".vtk");
                        } else {
                            append2VTK(result, variableName, outputPath + "results" + subName + ".vtk");
                        }
                        break;
                    default:
                        throw new UnsupportedOperationException("Format not supported");
                }
            }

            if (precision == 0) {
                Mat.save(PXY, outputPath + "vertices.txt");
            }

        } else {
            int listSize = names.size();
            for (int s = 0; s < listSize; s++) {
                ResultsCollection resCol = new ResultsCollection();
                String variableName = names.pop();
                System.out.println();
                System.out.println(variableName);
                System.out.println("|        |");
                double[] value = eqn.getResults(W[0], new double[eqn.nEqs() * dim], new double[dim], variableName);
                for (int i = 0; i < TP.length; i++) {
                    double[][] vertices = new double[TP[i].length][dim];
                    for (int j = 0; j < TP[i].length; j++) {
                        System.arraycopy(PXY[TP[i][j]], 0, vertices[j], 0, dim);
                    }

                    ElementType elemType = ElementType.elementTypeFactory(elemsType[i], order[i], par.volumeQuardatureOrder, par.faceQuardatureOrder);
                    Transformation transform = elemType.getVolumeTransformation(vertices, fCurv[i], par);

                    Basis basis = elemType.getBasis(transform);
                    int nBasis = basis.nBasis;
                    LocalElementSubdivision triLoc = new LocalElementSubdivision(elemsType[i], precision);
                    double[][] xCoords = transform.getX(triLoc.xiCoord);
                    double[][][] dXiBasis = basis.getDerBasisXi(triLoc.xiCoord, dim);
                    double[][][] dXBasis = transform.transformBasis(triLoc.xiCoord, dXiBasis);
                    double[][] values = new double[xCoords.length][value.length];
                    int nEqs = eqn.nEqs();
                    for (int j = 0; j < xCoords.length; j++) {
                        double[] Wpoint = new double[nEqs];
                        double[] dWpoint = new double[dim * nEqs];
                        for (int k = 0; k < Wpoint.length; k++) {
                            for (int m = 0; m < nBasis; m++) {
                                Wpoint[k] += W[i][k * nBasis + m] * basis.basisFun(m, triLoc.xiCoord[j]);
                                for (int d = 0; d < dim; d++) {
                                    dWpoint[nEqs * d + k] += W[i][k * nBasis + m] * dXBasis[j][m][d];
                                }
                            }
                        }
                        values[j] = eqn.getResults(Wpoint, dWpoint, xCoords[j], variableName);
                    }
                    resCol.add(i, xCoords, triLoc.TP, values, triLoc.localType);

                    if (TP.length > 10 && (i + 1) % (TP.length / 10) == 0) {
                        System.out.print("*");
                    }
                }
                resCol.saveResults(variableName, format, s);
            }
            System.out.println();
        }
    }

    private void save2VTK(double[][] result, double[][] PXY, int[][] TP, int[] elementType, String variableName, String fileName) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName))) {
            writer.write("# vtk DataFile Version 2.6");
            writer.newLine();
            writer.write("output file");
            writer.newLine();
            writer.write("ASCII");
            writer.newLine();
            writer.write("DATASET UNSTRUCTURED_GRID");
            writer.newLine();
            writer.newLine();

            writer.write("POINTS " + PXY.length + " float");
            writer.newLine();
            for (int i = 0; i < PXY.length; i++) {
                for (int j = 0; j < PXY[i].length; j++) {
                    writer.write(PXY[i][j] + " ");
                }
                for (int j = PXY[i].length + 1; j <= 3; j++) {
                    writer.write("0 ");
                }
                writer.newLine();
            }
            writer.newLine();

            int sumTP = 0;
            for (int i = 0; i < TP.length; i++) {
                for (int j = 0; j < TP[i].length; j++) {
                    sumTP++;
                }
            }
            writer.newLine();

            writer.write("CELLS " + TP.length + " " + (sumTP + TP.length));
            writer.newLine();
            for (int i = 0; i < TP.length; i++) {
                writer.write(TP[i].length + " ");
                for (int j = 0; j < TP[i].length; j++) {
                    writer.write(TP[i][j] + " ");
                }
                writer.newLine();
            }
            writer.newLine();

            writer.write("CELL_TYPES " + TP.length);
            writer.newLine();
            for (int i = 0; i < TP.length; i++) {
                writer.write(Integer.toString(VTKType(elementType[i])));
                writer.newLine();
            }
            writer.newLine();

            writer.write("POINT_DATA " + PXY.length);
            writer.newLine();
            if (result[0].length == 1) { // scalar data
                writer.write("SCALARS " + variableName + " float 1");
                writer.newLine();
                writer.write("LOOKUP_TABLE default");
                writer.newLine();
                for (int i = 0; i < PXY.length; i++) {
                    writer.write(Double.toString(result[i][0]));
                    writer.newLine();
                }
                writer.close();
            } else { // vector data
                writer.write("VECTORS " + variableName + " float");
                writer.newLine();
                for (int i = 0; i < PXY.length; i++) {
                    for (int k = 0; k < result[i].length; k++) {
                        writer.write(Double.toString(result[i][k]) + " ");
                    }
                    writer.newLine();
                }
                writer.newLine();
                writer.close();
            }
        }
    }

    private void append2VTK(double[][] result, String variableName, String fileName) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fileName, true))) {
            if (result[0].length == 1) { // scalar data
                writer.write("SCALARS " + variableName + " float 1");
                writer.newLine();
                writer.write("LOOKUP_TABLE default");
                writer.newLine();
                for (int i = 0; i < result.length; i++) {
                    writer.write(Double.toString(result[i][0]));
                    writer.newLine();
                }
                writer.close();
            } else { // vector data
                writer.write("VECTORS " + variableName + " float");
                writer.newLine();
                for (int i = 0; i < result.length; i++) {
                    for (int k = 0; k < result[i].length; k++) {
                        writer.write(Double.toString(result[i][k]) + " ");
                    }
                    writer.newLine();
                }
                writer.newLine();
                writer.close();
            }
        }
    }

    int VTKType(int elementType) {
        switch (firstDigit(elementType)) {
            case 1:
                return 1;
            case 2:
                return 3;
            case 3:
                return 5;
            case 4:
                return 9;
            case 5:
                return 10;
            case 6:
                return 12;
            case 7:
                return 13;
            case 8:
                return 14;
            default:
                return 0;
        }
    }

    private void smoothResults(double[][] r, int[][] TP, int smoothingIter) {
        int m = r.length;
        int n = r[0].length;
        for (int iter = 0; iter < smoothingIter; iter++) {
            double[][] rNew = new double[m][n];
            int[] sum = new int[m];
            for (int i = 0; i < TP.length; i++) {
                for (int j = 0; j < TP[i].length; j++) {
                    for (int k = 0; k < TP[i].length; k++) {
                        if (k != j) {
                            for (int var = 0; var < n; var++) {
                                rNew[TP[i][k]][var] += r[TP[i][j]][var];
                            }
                            sum[TP[i][k]] += 1;
                        }
                    }
                }
            }
            double alfa = 0.8;
            for (int i = 0; i < m; i++) {
                for (int var = 0; var < n; var++) {
                    r[i][var] = alfa * r[i][var] + (1 - alfa) * rNew[i][var] / sum[i];
                }
            }
        }
    }

    class LocalElementSubdivision {

        double[][] xiCoord;
        int[][] TP;
        int[] localType;

        LocalElementSubdivision(int elemType, int n) {
            int ne, s;
            double[] x;
            int[][] index;

            switch (firstDigit(elemType)) {
                case 2: // line
                    TP = new int[n][2];
                    localType = new int[n];
                    n = n + 1;
                    x = Mat.linspace(-1, 1, n);
                    xiCoord = new double[n][1];
                    for (int i = 0; i < n; i++) {
                        xiCoord[i][0] = x[i];
                    }
                    for (int i = 0; i < n - 1; i++) {
                        TP[i] = new int[]{i + 1, i};
                        localType[i] = 2;
                    }
                    break;

                case 3: // triangle
                    int nt = 1;
                    for (int i = 1; i < n; i++) {
                        nt = nt + 2 * i + 1;
                    }
                    TP = new int[nt][3];
                    localType = new int[nt];
                    n = n + 1;
                    index = new int[n][n];
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            index[i][j] = -1;
                        }
                    }
                    ne = n * (1 + n) / 2;
                    x = Mat.linspace(0, 1, n);
                    xiCoord = new double[ne][2];
                    s = 0;
                    for (int i = 0; i < n; i++) {
                        double[] y = Mat.linspace(0, x[i], i + 1);
                        for (int j = 0; j < i + 1; j++) {
                            xiCoord[s][0] = y[j];
                            xiCoord[s][1] = x[i] - y[j];
                            index[i][j] = s;
                            s = s + 1;
                        }
                    }

                    // structure generation
                    s = 0;
                    for (int i = 0; i < n - 1; i++) {
                        for (int j = 0; j < n - 1; j++) {
                            if (index[i][j] > -1 && index[i + 1][j] > -1 && index[i + 1][j + 1] > -1) {
                                TP[s] = new int[]{index[i][j], index[i + 1][j], index[i + 1][j + 1]};
                                localType[s] = 3;
                                s = s + 1;
                            }
                            if (index[i][j] > -1 && index[i + 1][j + 1] > -1 && index[i][j + 1] > -1) {
                                TP[s] = new int[]{index[i][j], index[i + 1][j + 1], index[i][j + 1]};
                                localType[s] = 3;
                                s = s + 1;
                            }
                        }
                    }
                    break;

                case 4: // quad
                    TP = new int[n * n][4];
                    localType = new int[n * n];
                    n = n + 1;
                    index = new int[n][n];
                    ne = n * n;
                    x = Mat.linspace(0, 1, n);
                    xiCoord = new double[ne][2];
                    s = 0;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            xiCoord[s][0] = x[i];
                            xiCoord[s][1] = x[j];
                            index[i][j] = s;
                            s = s + 1;
                        }
                    }

                    // structure generation
                    s = 0;
                    for (int i = 0; i < n - 1; i++) {
                        for (int j = 0; j < n - 1; j++) {
                            TP[s] = new int[]{index[i][j], index[i + 1][j], index[i + 1][j + 1], index[i][j + 1]};
                            localType[s] = 4;
                            s = s + 1;
                        }
                    }
                    break;

                case 6: // hexaheder
                    TP = new int[n * n * n][];
                    localType = new int[n * n * n];
                    n = n + 1;
                    int[][][] indexHexa = new int[n][n][n];
                    ne = n * n * n;
                    x = Mat.linspace(0, 1, n);
                    xiCoord = new double[ne][3];
                    s = 0;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            for (int k = 0; k < n; k++) {
                                xiCoord[s][0] = x[i];
                                xiCoord[s][1] = x[j];
                                xiCoord[s][2] = x[k];
                                indexHexa[i][j][k] = s;
                                s = s + 1;
                            }
                        }
                    }

                    // structure generation
                    s = 0;
                    for (int i = 0; i < n - 1; i++) {
                        for (int j = 0; j < n - 1; j++) {
                            for (int k = 0; k < n - 1; k++) {
                                TP[s] = new int[]{indexHexa[i][j][k], indexHexa[i + 1][j][k], indexHexa[i + 1][j + 1][k], indexHexa[i][j + 1][k], indexHexa[i][j][k + 1], indexHexa[i + 1][j][k + 1], indexHexa[i + 1][j + 1][k + 1], indexHexa[i][j + 1][k + 1]};
                                localType[s] = 6;
                                s = s + 1;
                            }
                        }
                    }
                    break;
            }
        }
    }

    class ResultsCollection {

        HashMap<Integer, double[][]> coordsMap;
        HashMap<Integer, int[][]> triLocMap;
        HashMap<Integer, double[][]> valuesMap;
        HashMap<Integer, int[]> localTypeMap;
        int s;
        int nCoord;
        int nTriLoc;
        int dimVal;
        int n;

        ResultsCollection() {
            coordsMap = new HashMap();
            triLocMap = new HashMap();
            valuesMap = new HashMap();
            localTypeMap = new HashMap();
            s = 0;
            nCoord = 0;
            nTriLoc = 0;
            dimVal = 0;
            n = 0;
        }

        void add(int i, double[][] xCoords, int[][] triLoc, double[][] values, int[] localType) {
            coordsMap.put(i, xCoords);
            triLoc = Mat.plus(triLoc, s);
            triLocMap.put(i, triLoc);
            valuesMap.put(i, values);
            localTypeMap.put(i, localType);
            s += xCoords.length;
            nCoord += xCoords.length;
            nTriLoc += triLoc.length;
            dimVal = values[0].length;
            if (dimVal > 1) {
                dimVal = 3; // correction for vector vtk format
            }
            n++;
        }

        // saving
        void saveResults(String variableName, String format, int append) throws IOException {
            double[][] vertices = new double[nCoord][dim];
            int p = 0;
            for (int i = 0; i < n; i++) {
                double[][] coords = coordsMap.get(i);
                for (int j = 0; j < coords.length; j++) {
                    for (int k = 0; k < coords[0].length; k++) {
                        vertices[p][k] = coords[j][k] * meshScale * lRef;
                    }
                    p++;
                }
            }

            int[][] tri = new int[nTriLoc][];
            p = 0;
            for (int i = 0; i < n; i++) {
                int[][] triLoc = triLocMap.get(i);
                for (int j = 0; j < triLoc.length; j++) {
                    tri[p] = new int[triLoc[0].length];
                    for (int k = 0; k < triLoc[0].length; k++) {
                        tri[p][k] = triLoc[j][k];
                    }
                    p++;
                }
            }

            double[][] val = new double[nCoord][dimVal];
            p = 0;
            for (int i = 0; i < n; i++) {
                double[][] value = valuesMap.get(i);
                for (int j = 0; j < value.length; j++) {
                    for (int k = 0; k < value[0].length; k++) {
                        val[p][k] = value[j][k];
                    }
                    p++;
                }
            }

            int[] localElementType = new int[nTriLoc];
            p = 0;
            for (int i = 0; i < n; i++) {
                int[] localType = localTypeMap.get(i);
                for (int j = 0; j < localType.length; j++) {
                    localElementType[p] = localType[j];
                    p++;
                }
            }

            // mesh saving
            File directory = new File(outputPath);
            if (!directory.exists()) {
                directory.mkdir();
            }
            switch (format) {
                case "txt":
                    Mat.save(vertices, outputPath + "vertices.txt");
                    Mat.save(tri, outputPath + "elements.txt");
                    Mat.save(val, outputPath + variableName + ".txt");
                    break;
                case "vtk":
                    if (append == 0) {
                        save2VTK(val, vertices, tri, localElementType, variableName, outputPath + "results.vtk");
                    } else {
                        append2VTK(val, variableName, outputPath + "results.vtk");
                    }
                    break;
                default:
                    throw new UnsupportedOperationException("Format not supported");
            }
        }
    }
}
