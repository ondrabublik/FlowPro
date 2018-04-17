package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import flowpro.core.Mesh.Element;
import flowpro.api.Functional;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class OptimisationToolExport {

    Mesh mesh;
    Functional fun;
    String simulationPath;
    String optimisationPath;
    String mode;
    public static final String PARAMETER_FILE_NAME = "parameters.txt";

    OptimisationToolExport(Solver solver, String simulationPath, String mode) throws IOException {
        this.mesh = solver.getMesh();
        this.simulationPath = simulationPath;
        this.mode = mode;
        optimisationPath = simulationPath + "/optimisation/";
    }

    public void export() throws IOException {
        FlowProProperties props = new FlowProProperties();
        props.load(new FileInputStream(simulationPath + PARAMETER_FILE_NAME));
        System.out.println("Functional of interest:" + props.getString("functionalOfInterest"));
        int numberOfAlpha = props.getInt("numberOfAlpha");
        System.out.println("Number of optimalization parameters: " + numberOfAlpha);

        fun = (new FunctionalFactory()).getFunctional(simulationPath + PARAMETER_FILE_NAME);
        for (Element elem : mesh.getElems()) {
            elem.optimalisationFunctional = fun;
        }

        double meshScale = 1;
        if (props.containsKey("meshScale")) {
            meshScale = props.getDouble("meshScale");
            System.out.println("Mesh scale parameter: " + meshScale);
        }

        // exporting
        File directory = new File(optimisationPath);
        if (!directory.exists()) {
            directory.mkdir();
        }

        exportFunctional(0);

        if (mode.equalsIgnoreCase("all")) {
            exportJacobiTranspose();
            exportMassMatrix();
            exportResiduum(0);
            exportFunctionalDerivative();

            for (int alpha = 0; alpha < numberOfAlpha; alpha++) {
                try {
                    Element[] elems = mesh.getElems();
                    double[][] PXY = Mat.loadDoubleMatrix(optimisationPath + "vertices" + (alpha + 1) + ".txt"); // mesh vertices coordinates
                    if (meshScale != 1) {
                        for (int i = 0; i < PXY.length; i++) {
                            for (int j = 0; j < PXY[i].length; j++) {
                                PXY[i][j] *= meshScale;
                            }
                        }
                    }
                    for (int i = 0; i < mesh.nElems; i++) {
                        for (int j = 0; j < elems[i].vertices.length; j++) {
                            System.arraycopy(PXY[elems[i].TP[j]], 0, elems[i].vertices[j], 0, mesh.getEqn().dim());
                        }
                    }
                    mesh.init();
                    exportFunctional(alpha + 1);
                    exportResiduum(alpha + 1);
                } catch (Exception e) {
                    System.out.println("Error, file " + optimisationPath + "vertices" + (alpha + 1) + ".txt not found!");
                }
            }

            // export sparse vector print
            exportPrintMatrix();
        }
    }

    public void exportJacobiTranspose() throws IOException {
        Element[] elems = mesh.getElems();
        FileWriter fw = new FileWriter(optimisationPath + "J.txt");
        try (BufferedWriter out = new BufferedWriter(fw)) {
            for (int k = 0; k < mesh.nElems; k++) {
                elems[k].exportLocalJacobiMatrix();
                int n = elems[k].nBasis * mesh.nEqs;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        out.write(elems[k].gi_U[i] + " " + elems[k].gi_U[j] + " " + elems[k].ADiag[j][i]);
                        out.newLine();
                    }
                }
                for (int face = 0; face < elems[k].nFaces; face++) {
                    if (elems[k].TT[face] > -1) {
                        Element elemR = elems[elems[k].TT[face]];
                        for (int i = 0; i < n; i++) {
                            for (int j = 0; j < elemR.nBasis * mesh.nEqs; j++) {
                                out.write(elems[k].gi_U[i] + " " + elemR.gi_U[j] + " " + elems[k].ANeighs[face].MR[j][i]);
                                out.newLine();
                            }
                        }
                    }
                }
            }
            out.close();
        }
    }

    public void exportResiduum(int alpha) throws IOException {
        Element[] elems = mesh.getElems();
        FileWriter fw = new FileWriter(optimisationPath + "R" + alpha + ".txt");
        try (BufferedWriter out = new BufferedWriter(fw)) {
            for (int i = 0; i < mesh.nElems; i++) {
                elems[i].exportLocalR();
                for (int j = 0; j < elems[i].RHS_loc.length; j++) {
                    out.write(elems[i].RHS_loc[j] + " ");
                    out.newLine();
                }
            }
            out.close();
        }
    }

    public void exportFunctionalDerivative() throws IOException {
        Element[] elems = mesh.getElems();
        FileWriter fw = new FileWriter(optimisationPath + "dIdw.txt");
        try (BufferedWriter out = new BufferedWriter(fw)) {
            for (int i = 0; i < mesh.nElems; i++) {
                elems[i].exportLocalFunctionalDerivative();
                for (int j = 0; j < elems[i].optimFunDer.length; j++) {
                    out.write(elems[i].optimFunDer[j] + " ");
                    out.newLine();
                }
            }
            out.close();
        }
    }

    public void exportFunctional(int alpha) throws IOException {
        Element[] elems = mesh.getElems();
        double f = 0;
        for (Element elem : elems) {
            f += elem.computeFunctional(null);
        }
        FileWriter fw = new FileWriter(optimisationPath + "I" + alpha + ".txt");
        try (BufferedWriter out = new BufferedWriter(fw)) {
            out.write(f + " ");
            out.newLine();
            out.close();
        }
    }

    public void exportMassMatrix() throws IOException {
        Element[] elems = mesh.getElems();
        FileWriter fw = new FileWriter(optimisationPath + "M.txt");
        try (BufferedWriter out = new BufferedWriter(fw)) {
            for (int k = 0; k < mesh.nElems; k++) {
                int nb = elems[k].nBasis;
                int nEq = mesh.nEqs;
                for (int m = 0; m < nEq; m++) {
                    for (int i = 0; i < nb; i++) {
                        for (int j = 0; j < nb; j++) {
                            out.write(elems[k].gi_U[nb * m + i] + " " + elems[k].gi_U[nb * m + j] + " " + elems[k].M[i][j]);
                            out.newLine();
                        }
                    }
                }
            }
            out.close();
        }
    }

    void exportPrintMatrix() throws IOException {
        Element[] elems = mesh.getElems();
        FileWriter fw = new FileWriter(optimisationPath + "printMatrix.txt");
        try (BufferedWriter out = new BufferedWriter(fw)) {
            for (int k = 0; k < mesh.nElems; k++) {
                for (int i = 0; i < elems[k].nBasis; i++) {
                    for (int j = 0; j < mesh.nEqs; j++) {
                        out.write((k * mesh.nEqs + j) + " " + elems[k].gi_U[j * elems[k].nBasis + i] + " " + elems[k].Is[i]);
                        out.newLine();
                    }
                }
            }
            out.close();
        }
    }
}
