package flowpro.core.meshDeformation;

import flowpro.api.*;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import flowpro.core.elementType.*;
import java.io.*;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author obublik
 */
public abstract class Deformation implements Serializable {

    private static final Logger LOG = LoggerFactory.getLogger(Deformation.class);
    public int nBodies;
    public double[][] center;
    public double[][] blendFuns;
    final public Equation eqn;
    public int nUserValues = 2;
    public double[][] userDefTerm; // user deffined term (important for some types of application)

    public Deformation(Parameters par, Equation eqn, int[][] TEale) throws IOException {
        this.eqn = eqn;
        nBodies = countBodies(TEale);

        if (nBodies > 0) {
            LOG.info("Number of bodies: " + nBodies);
        } else {
            LOG.warn("Moving structure not found! movingMesh is set to false!");
            //par.movingMesh = false;
        }
    }

    abstract public void newMeshPositionAndVelocity(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov);

    abstract public void nextTimeLevel(Element[] elems);

    abstract public void calculateForces(Element[] elems, MeshMove[] mshMov);

    public void recalculateMesh(Element[] elems) { // dodelat paralelni verzi
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.transform.computeTransform(elem.vertices);
            }
        }
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.Int.recalculateVolumeIntegrals(elem);
            }
        }
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                elem.Int.recalculateFaceIntegrals(elem);
                elem.computeGeometry();
                elem.recalculateMassMatrix();
            }
        }
    }

    public void relaxFirstIteration(Element[] elems, double dt) { // prevent solution error, when jump in mesh position in first iteration occure
        for (Element elem : elems) {
            if (elem.insideComputeDomain) {
                for (int j = 0; j < elem.nVertices; j++) {
                    for (int d = 0; d < elem.dim; d++) {
                        elem.verticesOld2[j][d] = elem.vertices[j][d] - 2*dt*elem.Uinit[j][d];
                        elem.verticesOld[j][d] = elem.vertices[j][d] - dt*elem.Uinit[j][d];
                        elem.U[j][d] = elem.Uinit[j][d];
                    }
                }
            }
        }
    }

    public double[][] getBlendFuns() {
        return blendFuns;
    }

    private int countBodies(int[][] TEale) {
        nBodies = 0;
        for (int i = 0; i < TEale.length; i++) {
            for (int j = 0; j < TEale[i].length; j++) {
                if (TEale[i][j] - 1 > nBodies) {
                    nBodies = TEale[i][j] - 1;
                }
            }
        }
        return nBodies;
    }

    abstract public FluidForces[] getFluidForces();

    public double[][] getUserDef() {
        return userDefTerm;
    }

    public void setCenters(double[][] center) {
        this.center = center;
    }

    public void calculateBlendingFunctions(double[][] PXY, int[][] TP, int[][] TT, int[][] TEale, int[] elemsType, String meshPath) throws IOException {

        try {
            blendFuns = Mat.loadDoubleMatrix(meshPath + "blendingFunctions.txt");
        } catch (FileNotFoundException ex) {
            System.out.println("Generating blending functions!");
            ElementType[] eTyp = new ElementType[8];
            eTyp[1] = new pointElement(1,1,1);
            eTyp[2] = new lineElement(1,1,1);
            eTyp[3] = new triangleElement(1,1,1);
            eTyp[4] = new squareElement(1,1,1);
            eTyp[5] = new tetrahedralElement(1,1,1);
            eTyp[6] = new hexahedralElement(1,1,1);
            eTyp[7] = new prismElement(1,1,1);

            int nPoints = PXY.length;
            int nElem = TP.length;
            int dim = PXY[0].length;
            blendFuns = new double[nPoints][nBodies];

            // calculating sparse matrix size
            int nA = 0;
            for (int i = 0; i < nElem; i++) {
                for (int j = 0; j < eTyp[elemsType[i]].nFaces; j++) {
                    if (dim == 2) {
                        nA = nA + 6;
                    } else {
                        nA = nA + 4 * eTyp[elemsType[i]].getFaceIndexes(j).length;
                    }
                }
            }

            int[] IA = new int[nA];
            int[] JA = new int[nA];
            double[] HA = new double[nA];
            double[] sumHA = new double[nPoints];
            int s = 0;
            for (int i = 0; i < nElem; i++) {
                for (int j = 0; j < eTyp[elemsType[i]].nFaces; j++) {
                    int[] faceIndexes = eTyp[elemsType[i]].getFaceIndexes(j);
                    int nIndexes = faceIndexes.length;
                    int nk = nIndexes;
                    if (dim == 2) {
                        nk = nIndexes - 1;
                    }
                    double koef = 1;
                    if (TEale[i][j] != 0) {
                        koef = 2;
                    }
                    for (int k = 0; k < nk; k++) {
                        int kp = (k + 1) % nIndexes;
                        double S = 0;
                        for (int d = 0; d < dim; d++) {
                            S += (PXY[TP[i][faceIndexes[kp]]][d] - PXY[TP[i][faceIndexes[k]]][d]) * (PXY[TP[i][faceIndexes[kp]]][d] - PXY[TP[i][faceIndexes[k]]][d]);
                        }
                        S = Math.sqrt(S);

                        IA[s] = TP[i][faceIndexes[k]];
                        JA[s] = TP[i][faceIndexes[k]];
                        HA[s] = -koef / S;
                        s++;

                        IA[s] = TP[i][faceIndexes[k]];
                        JA[s] = TP[i][faceIndexes[kp]];
                        HA[s] = koef / S;
                        s++;

                        IA[s] = TP[i][faceIndexes[kp]];
                        JA[s] = TP[i][faceIndexes[kp]];
                        HA[s] = -koef / S;
                        s++;

                        IA[s] = TP[i][faceIndexes[kp]];
                        JA[s] = TP[i][faceIndexes[k]];
                        HA[s] = koef / S;
                        s++;

                        sumHA[TP[i][faceIndexes[k]]] += koef / S;
                        sumHA[TP[i][faceIndexes[kp]]] += koef / S;
                    }
                }
            }
            for (int i = 0; i < IA.length; i++) {
                HA[i] /= sumHA[IA[i]];
            }

            // calculating
            boolean[] boundary = new boolean[nPoints];
            Arrays.fill(boundary, false);
            for (int i = 0; i < nElem; i++) {
                for (int j = 0; j < eTyp[elemsType[i]].nFaces; j++) {
                    int[] faceIndexes = eTyp[elemsType[i]].getFaceIndexes(j);
                    if (TEale[i][j] > 0) {
                        for (int k = 0; k < faceIndexes.length; k++) {
                            boundary[TP[i][faceIndexes[k]]] = true;
                        }
                    }
                }
            }
            for (int body = 0; body < nBodies; body++) {
                double[] x = new double[nPoints];
                for (int i = 0; i < nElem; i++) {
                    for (int j = 0; j < eTyp[elemsType[i]].nFaces; j++) {
                        int[] faceIndexes = eTyp[elemsType[i]].getFaceIndexes(j);
                        if (TEale[i][j] == body + 2) {
                            for (int k = 0; k < faceIndexes.length; k++) {
                                x[TP[i][faceIndexes[k]]] = 1;
                            }
                        }
                    }
                }
                calculate(IA, JA, HA, x, boundary, body);
            }

            Mat.save(blendFuns, meshPath + "blendingFunctions.txt");
        }
    }

    public void calculate(int[] IA, int[] JA, double[] HA, double[] x, boolean[] boundary, int body) {
        // Gauss-Seidel method
        double[] xNew = new double[x.length];

        for (int iter = 0; iter < 10000; iter++) {
            for (int i = 0; i < IA.length; i++) {
                xNew[IA[i]] += HA[i] * x[JA[i]];
            }
            double res = 0;
            for (int i = 0; i < x.length; i++) {
                if (!boundary[i]) {
                    if (Math.abs(xNew[i] - x[i]) > res) {
                        res = Math.abs(xNew[i] - x[i]);
                    }
                    x[i] = xNew[i];
                }
            }
            if (res < 1e-6) {
                System.out.println("Blending function of " + (body + 1) + "-th body converged in " + iter + " iteration.");
                break;
            }
        }

        for (int i = 0; i < x.length; i++) {
            blendFuns[i][body] = x[i];
        }
    }
}
