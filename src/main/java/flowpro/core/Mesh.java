package flowpro.core;

import flowpro.core.meshDeformation.*;
import flowpro.core.quadrature.QuadratureCentral;
import flowpro.api.Equation;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.elementType.ElementType;
import flowpro.api.SolutionMonitor;
import flowpro.core.parallel.Domain.Subdomain;
import flowpro.core.element.*;
import flowpro.core.solver.MasterSolver;
import java.io.*;
import java.util.Arrays;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * implementace DGFEM
 */
public class Mesh implements Serializable {

    private static final Logger LOG = LoggerFactory.getLogger(Mesh.class);

    private final Equation eqn;
    private final Deformation dfm;
    public final int nEqs;        // number of equations
    private final Parameters par;  // parameters for the solver
    QuadratureCentral qRules;     // set of quadrature rules
    private final Element[] elems;   // array of all elements in the mesh
    public final int nElems;
    public final int nPoints;
    public int dofs;   // number of degrees of freedom

    /* atributes for the prallel mode */
    private int[] boundary;
    public int[] save;
    public double t; // time

    // domain monitor
    public SolutionMonitor solMonitor;
    public double[] integralMonitor;

    /**
     * Creates an instance of class Element for each element in the mesh and
     * places them into an array (serves as Element factory). Method getElems()
     * returns the array of elements. This constructor is intended to be used
     * only in the parallel mode.
     *
     * @param eqn
     * @param dfm
     * @param par
     * @param solMonitor
     * @param qRules
     * @param PXY
     * @param UXY
     * @param elemsOrder
     * @param wallDistance
     * @param externalField
     * @param elemsType
     * @param TP
     * @param TT
     * @param TEale
     * @param TEshift
     * @param shift
     * @param initW
     * @param fCurv
     * @param domain
     * @throws IOException
     */
    public Mesh(Equation eqn, Deformation dfm, Parameters par, SolutionMonitor solMonitor, QuadratureCentral qRules, double[][] PXY, double[][] UXY, int[] elemsOrder, double[] wallDistance, double[][] externalField, int[] elemsType, int[][] TP, int[][] TT, int[][] TEale, int[][] TEshift, double[][] shift,
            FaceCurvature[] fCurv, double[][] initW, Subdomain domain) throws IOException {
        this.eqn = eqn;
        this.dfm = dfm;
        this.nEqs = eqn.nEqs();
        this.par = par;
        this.solMonitor = solMonitor;
        this.qRules = qRules;
        this.boundary = domain.boundary;
        this.save = domain.save;

        // vytvoreni vypocetnich elementu
        nElems = TT.length;
        nPoints = PXY.length;
        elems = new Element[nElems];
        for (int i = 0; i < nElems; i++) {
            int nVertices = TP[i].length;
            int nEdges = TT[i].length;
            double[][] vertices = new double[nVertices][eqn.dim()];
            double[][] meshVelocity = new double[nVertices][eqn.dim()];
            double[] wallDistancee = new double[nVertices];
            int[] TTe = new int[nEdges];
            int[] TPe = new int[nVertices];
            int[] TEalee = new int[nEdges];
            int[] TEshifte = new int[nEdges * eqn.dim()];
            for (int j = 0; j < nVertices; j++) {
                System.arraycopy(PXY[TP[i][j]], 0, vertices[j], 0, eqn.dim());
                wallDistancee[j] = wallDistance[TP[i][j]];
                TPe[j] = TP[i][j];
            }
            if (UXY != null) {
                for (int j = 0; j < nVertices; j++) {
                    System.arraycopy(UXY[TP[i][j]], 0, meshVelocity[j], 0, eqn.dim());
                }
            }
            for (int j = 0; j < nEdges; j++) {
                TTe[j] = TT[i][j];
                TEalee[j] = TEale[i][j];
                TEshifte[j] = TEshift[i][j];
            }

            double[][] externalFielde = null;
            if (par.externalField) {
                int nExternalField = externalField[0].length;
                externalFielde = new double[nVertices][nExternalField];
                for (int j = 0; j < nVertices; j++) {
                    System.arraycopy(externalField[TP[i][j]], 0, externalFielde[j], 0, nExternalField);
                }
            }

            // read or compute blending function
            double[][] blendFunse = null;
            if (par.movingMesh) {
                double[][] blendFuns = dfm.getBlendFuns();
                blendFunse = new double[nVertices][dfm.nBodies];
                for (int j = 0; j < nVertices; j++) {
                    for (int k = 0; k < dfm.nBodies; k++) {
                        blendFunse[j][k] = blendFuns[TP[i][j]][k];
                    }
                }
            }

            // create element type
            elems[i] = getSpatialMethod(par.spatialMethod);

            // seting maximal allowed spatial order for choosen method
            int elemOrder = Math.min(elemsOrder[i], elems[i].maxMethodSpatialOrder);
            int volumeQuardatureOrder = Math.min(par.volumeQuardatureOrder, elems[i].maxMethodSpatialOrder);
            int faceQuardatureOrder = Math.min(par.faceQuardatureOrder, elems[i].maxMethodSpatialOrder);
            par.order = Math.min(par.order, elems[i].maxMethodSpatialOrder);

            if (i == 0 && elemOrder != elemsOrder[i]) {
                LOG.warn("method does not support spatial order higher than " + elemOrder);
            }

            ElementType elemType = ElementType.elementTypeFactory(elemsType[i], elemOrder, volumeQuardatureOrder, faceQuardatureOrder);
            elems[i].set(i, vertices, meshVelocity, wallDistancee, externalFielde, TTe, TPe, TEalee, TEshifte, shift, fCurv[i], blendFunse, initW[i], this, elemType);
        }

        // seting elements position in domain (load and insideComputeDomain parts)
        for (int i : boundary) {
            elems[i].insideComputeDomain = false;
        }
        for (int i : domain.interior) {
            elems[i].insideMetisDomain = true;
        }
        for (int i : domain.save1) {
            elems[i].gmresSave = true;
        }
    }

    public void init() throws IOException {
        // Firstly, init basis function at each of elements
        for (Element elem : elems) {
            elem.initBasis();
        }
        System.out.println("basis       **********");
        System.out.print("integration ");

        // Secondly init volume integration rules on each of elements
        for (int i = 0; i < nElems; i++) {
            elems[i].initIntegration();
            if (nElems > 10 && i % (nElems / 10) == 0) {
                System.out.print("*");
            }
        }
        System.out.println();

        System.out.print("geometry    ");
        int dofs0 = 0;
        int nBrokenElements = 0;
        for (int i = 0; i < nElems; i++) {
            dofs0 = elems[i].setGlobalIndex(dofs0);
            if (elems[i].insideComputeDomain) {
                elems[i].Int.initNeighbours(elems, elems[i]);   // init face integration rules
                elems[i].computeGeometry();                     // compute geometrical relations
                elems[i].computeMassMatrixAndBasisWeights();                 // compute mass matrix
                elems[i].initMethod(par.props);        // for damping
                elems[i].createTimeIntegrationElementObject(elems[i]);

                // checking geometry
                boolean isOK = elems[i].geometryCheck(false);
                if (!isOK) {
                    nBrokenElements++;
                }
            }
            elems[i].insideAssemblerDomain = elems[i].insideMetisDomain;
            for (int j = 0; j < elems[i].nFaces; j++) {
                if (elems[i].TT[j] > -1 && !elems[elems[i].TT[j]].insideMetisDomain) {
                    elems[i].insideAssemblerDomain = false;
                }
            }
            if (nElems > 10 && i % (nElems / 10) == 0) {
                System.out.print("*");
            }
        }
        dofs = dofs0;

        System.out.println();
        LOG.info("Degrese of freedom: " + dofs);
        // initial condition on each of elements
        for (Element elem : elems) {
            elem.initCondition();
        }

        if (nBrokenElements > 0) {
            LOG.error("Broken elements = " + nBrokenElements);
        } else {
            LOG.info("Broken elements = 0 ");
        }
    }

    /**
     *
     * @return array of all elements in the mesh
     */
    public Element[] getElems() {
        return elems;
    }

    public Parameters getPar() {
        return par;
    }

    public Equation getEqn() {
        return eqn;
    }

    public Deformation getDfm() {
        return dfm;
    }

    public QuadratureCentral getQRules() {
        return qRules;
    }

    public Solution getSolution() {
        if (par.movingMesh) {
            return new Solution(getW(), getAvgW(), getMeshPosition(), getMeshVelocity());
        } else {
            return new Solution(getW(), getAvgW(), getMeshPosition());
        }
    }

    /**
     *
     * @return nElems x 4 array containing integral averages of the solution on
     * all elements
     */
    public double[][] getAvgW() {
        double[][] avgW = new double[nElems][];
        for (int i = 0; i < nElems; ++i) {
            if (elems[i].insideComputeDomain) {
                avgW[i] = elems[i].calculateAvgW();
            }
        }
        return avgW;
    }

    /**
     *
     * @return 2D array whose rows are basis coefficients of the corresponding
     * element
     */
    public double[][] getW() {
        double[][] W = new double[nElems][];
        for (int i = 0; i < nElems; ++i) {
            Element elem = elems[i];
            W[i] = new double[nEqs * elem.nBasis];
            for (int j = 0; j < elem.nBasis; ++j) {
                for (int k = 0; k < nEqs; ++k) {
                    W[i][k * elem.nBasis + j] = elem.W[k * elem.nBasis + j];
                }
            }
        }
        return W;
    }

    /**
     *
     * return array containing integral averages of the solution through all
     * elements
     *
     */
    public void computeSolutionMonitor() {
        int nIntegralMonitor = solMonitor.getNumberOfMonitoredValues();
        integralMonitor = new double[nIntegralMonitor];
        for (Element elem : elems) {
            if (elem.insideMetisDomain) {
                double[] aux = elem.computeIntegralSolutionMonitor(nIntegralMonitor);
                for (int i = 0; i < nIntegralMonitor; i++) {
                    integralMonitor[i] += aux[i];
                }
            }
        }
        // combination of computed values
        solMonitor.combineMonitoredValues(integralMonitor);
    }

    public double[][] getMeshPosition() {
        double[][] PXY = new double[nPoints][eqn.dim()];
        for (Element elem : elems) {
            for (int j = 0; j < elem.TP.length; j++) {
                for (int d = 0; d < eqn.dim(); d++) {
                    PXY[elem.TP[j]][d] = elem.vertices[j][d];
                }
            }
        }
        return PXY;
    }

    public double[][] getMeshVelocity() {
        double[][] U = new double[nPoints][eqn.dim()];
        for (Element elem : elems) {
            for (int j = 0; j < elem.TP.length; j++) {
                for (int d = 0; d < eqn.dim(); d++) {
                    U[elem.TP[j]][d] = elem.U[j][d];
                }
            }
        }
        return U;
    }

    public double[][] getArtificialViscosity() {
        double[][] artVis = new double[nElems][nEqs];
        for (int i = 0; i < nElems; ++i) {
            for (int m = 0; m < nEqs; m++) {
                artVis[i][m] = elems[i].eps + elems[i].dampInner[m];
            }
        }
        return artVis;
    }

    public void updateTime(double dt) {
        t = t + dt;
    }

    public enum SpatialMethodType {
        DG, DGjacobi, DGpure, DGpureFast, FVM, DGpureIncompressible;

        public static void help() {
            System.out.println("********************************");
            System.out.println("HELP for parameter spatialMethod");
            System.out.println("list of possible values:");
            System.out.println(Arrays.asList(SpatialMethodType.values()));
            System.out.println("********************************");
        }
    }

    public Element getSpatialMethod(String spatialMethodName) throws IOException {
        Element elem = null;
        try {
            SpatialMethodType spatialMethod = SpatialMethodType.valueOf(spatialMethodName);
            switch (spatialMethod) {
                case DG:
                    elem = new DG();
                    break;
                case DGjacobi:
                    elem = new DGjacobi();
                    break;
                case DGpure:
                    elem = new DGpure();
                    break;
                case DGpureFast:
                    elem = new DGpureFast();
                    break;
                case DGpureIncompressible:
                    elem = new DGpureIncompressible();
                    break;
                case FVM:
                    elem = new FVM();
                    break;
            }
        } catch (IllegalArgumentException ex) {
            SpatialMethodType.help();
            throw new IOException("unknown spatial method " + spatialMethodName);
        }
        return elem;
    }
}
