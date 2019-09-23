/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.ElementData;
import flowpro.api.Equation;
import flowpro.api.Functional;
import flowpro.api.Mat;
import flowpro.core.Integration;
import flowpro.core.Mesh;
import flowpro.core.Parameters;
import flowpro.core.basis.Basis;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.elementType.ElementType;
import flowpro.core.quadrature.QuadratureCentral;
import flowpro.core.transformation.Transformation;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
/**
 * Represents an element in the mesh.
 */
public abstract class Element implements Serializable {

    // mesh
    private final Mesh mesh;
    // parameters
    Parameters par;
    // equations
    int nEqs;
    Equation eqn;
    // mesh elements
    Element[] elems;

    // time integration
    public TimeIntegrationElement ti;

    /* constants */
    public double c_IP;  // ???
    public double eps;  // ???

    /* geometry */
    public final int[] TT;     // indexes of neighbours
    public final int[] TP;
    public final int[] TEale;     // type of boundary for ALE
    public final int[] TEshift;
    public double[][] shift;
    public FaceCurvature fCurv;
    public int nFaces;       // number of faces
    public final int nVertices;       // number of vertices
    public double[][] vertices;  // components of vertices [xi,yi,zi]
    public double[][] verticesOld;
    public double[][] verticesOld2;
    public double[][] vertices0;  // components of vertices [xi,yi,zi]
    public double[][] U;  // x-components of velocity [ui, vi, wi]
    public double[][] Uinit;
    public double[][][] n;  // x-components of normals
    public double[] S;
    public double area;  // element area
    public double[] wallDistance; // vzdalennost od steny
    public double[] Xs;  // coordinates of element's midpoint
    public double[][] Xes; // coordinates of element's edge midpoint
    // size of element defined as 2 times the diameter of inscribed cyrcle
    public double elemSize;
    public boolean insideComputeDomain; // true if element is insideComputeDomain subdomain, false on the boundary
    public boolean insideMetisDomain; // true insideComputeDomain computational subdomain (producet by metis)
    public boolean insideAssemblerDomain; // interior for Jacobi matrix assemble
    public boolean gmresLoad, gmresSave;
    public double[][] blendFun; // blending function used for calculation of new mesh position

    // mesh
    public final int index;
    public int[] faceIndexReverse;

    public final int dim;
    public int nBasis; 	  // pocet bazovych funkci
    public Transformation transform;
    public Basis basis;
    public Integration Int;
    QuadratureCentral qRules;
    public ElementType elemType;

    //private double k_el; // kvalita elementu
    public double[] Is; // integral bazove funkce
    public double[][] M;	 // matice hmotnosti
    double[][] iM; // inverze matice hmotnosti (pouze pro implicitni metodu)
    double[][] Mo;	 // matice hmotnosti v predchozi casove hladine
    double[][] Mo2;	 // matice hmotnosti

    // damping
    double[][] TrunOrd;
    public double[] dampInner;
    public double[] innerIndicatorOld;

    // parametery pro prenos do equations
    public double[] centreVolumeInterpolant;
    public final ElementData elemData;

    // optimalisation
    public Functional optimalisationFunctional; // functional for optimalisation
    double[] optimFunDer;

    // matice globalnich indexu a globalni matice s pravou stranou
    public int[] gi_U;

    double[] initW;
    public double[] W;     // hodnota W v n+1 te casove hladine
    public double[] Wo;    // hodnota v n te casove hladine
    public double[] Wo2;   // hodnota v n-1 casove hladine

    // derivative correction
    double[] R;

    // external field
    double[][] externalField;

    public Element(int index, double[][] vertices, double[][] Uinit, double[] wallDistance, double[][] externalField, int[] TT, int[] TP, int[] TEale, int[] TEshift, double[][] shift, FaceCurvature fCurv, double[][] blendFun, double[] initW,
            Mesh mesh, ElementType elemType) throws IOException {

        this.mesh = mesh;
        eqn = mesh.getEqn();
        dim = eqn.dim();
        nEqs = eqn.nEqs();
        par = mesh.getPar();
        qRules = mesh.getQRules();
        elems = mesh.getElems();

        this.index = index;
        this.TT = TT;
        this.TP = TP;
        this.TEale = TEale;
        this.TEshift = TEshift;
        this.shift = shift;
        this.fCurv = fCurv;
        this.vertices = vertices;
        this.Uinit = Uinit;
        this.blendFun = blendFun;
        this.wallDistance = wallDistance;
        this.externalField = externalField;
        this.nFaces = elemType.nFaces;
        this.initW = initW;
        this.nVertices = vertices.length;
        this.elemType = elemType;

        insideComputeDomain = true;
        gmresLoad = false;
        gmresSave = false;

        U = new double[nVertices][dim];
        verticesOld = new double[nVertices][dim];
        verticesOld2 = new double[nVertices][dim];

        elemData = new ElementData(dim);
        centreVolumeInterpolant = new double[nVertices];
        Arrays.fill(centreVolumeInterpolant, 1.0 / nVertices);

        // damping
        dampInner = new double[nEqs];
        innerIndicatorOld = new double[nEqs];
        
        vertices0 = new double[nVertices][dim];
        for (int d = 0; d < dim; d++) {
            for (int i = 0; i < nVertices; i++) {
                vertices0[i][d] = vertices[i][d];
                U[i][d] = Uinit[i][d];
            }
        }
    }

    abstract public void initMethod();

    abstract public void initCondition();

    abstract public void computeMassMatrixAndBasisWeights();

    abstract public void recalculateMassMatrix();

    abstract public void residuum(double[] V, double[] K, double[][] KR);

    abstract public void residuumWall(int k, double[] V, double[] K, double[] KR);

    abstract public void limiter(boolean isFirstIter);

    public void createTimeIntegration(Element elem) throws IOException {
        ti = getTimeIntegrationElement(par.timeMethod, this);
        ti.set(elem);
        ti.init();
    }

    public void initBasis() throws IOException {
        transform = elemType.getVolumeTransformation(vertices, fCurv, par);
        basis = elemType.getBasis(transform);
        nBasis = basis.nBasis;
        Mo = new double[nBasis][nBasis];
        Mo2 = new double[nBasis][nBasis];
    }

    public void initIntegration() throws IOException {
        Int = new Integration(elemType, dim, basis, transform, TT, TEshift, shift, qRules);
    }

    public int nastav_globalni_index_U(int s) {
        gi_U = new int[nEqs * nBasis];
        for (int i = 0; i < nBasis * nEqs; i++) {
            gi_U[i] = s;
            s = s + 1;
        }
        return s;
    }

    public void computeInvertMassMatrix() {
        // vypocitava inverzi diagonaly, ktera se pouzije pro predpodminovac
        double[][] invM = Mat.invert(M);
        iM = new double[nEqs * nBasis][nEqs * nBasis];
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    iM[m * nBasis + i][m * nBasis + j] = invM[i][j];
                }
            }
        }
    }

    double[] interpolateVelocityAndFillElementDataObjectOnVolume(double[] innerInterpolant) {
        // interpolation of mesh velocity, and other data
        double[] u = new double[dim];
        double[] currentXi = new double[dim];
        elemData.currentWallDistance = 0;
        for (int j = 0; j < nVertices; j++) {
            for (int d = 0; d < dim; d++) {
                u[d] += innerInterpolant[j] * U[j][d];
                //elemData.currentX[d] += innerInterpolant[j] * vertices[j][d];
                currentXi[d] += innerInterpolant[j] * transform.coordsXi[j][d];
            }
            elemData.currentWallDistance += innerInterpolant[j] * wallDistance[j];
        }
        elemData.currentX = transform.getX(currentXi);

        if (par.externalField) {
            int nExternalField = externalField[0].length;
            elemData.externalField = new double[nExternalField];
            for (int j = 0; j < nVertices; j++) {
                for (int r = 0; r < nExternalField; r++) {
                    elemData.externalField[r] += innerInterpolant[j] * externalField[j][r];
                }
            }
        }
        elemData.currentT = mesh.t;
        elemData.integralMonitor = mesh.integralMonitor;
        elemData.elemIndex = index;
        if (par.solutionAverage) {
            elemData.Wavg = calculateAvgW();
        }

        return u;
    }

    double[] interpolateVelocityAndFillElementDataObjectOnFace(int k, double[] innerInterpolant, int[] edgeIndex) {
        // interpolation of mesh velocity
        double[] u = new double[dim];
        double[] currentXi = new double[dim];
        elemData.currentWallDistance = 0;
        for (int j = 0; j < Int.faces[k].nVerticesEdge; j++) {
            for (int d = 0; d < dim; d++) {
                u[d] += innerInterpolant[j] * U[edgeIndex[j]][d];
                elemData.meshVelocity[d] = u[d];
                //elemData.currentX[d] += innerInterpolant[j] * vertices[edgeIndex[j]][d];
                currentXi[d] += innerInterpolant[j] * transform.coordsXi[edgeIndex[j]][d];
            }
            elemData.currentWallDistance += innerInterpolant[j] * wallDistance[edgeIndex[j]];
        }
        elemData.currentX = transform.getX(currentXi);

        if (par.externalField) {
            int nExternalField = externalField[0].length;
            elemData.externalField = new double[nExternalField];
            for (int j = 0; j < Int.faces[k].nVerticesEdge; j++) {
                for (int r = 0; r < nExternalField; r++) {
                    elemData.externalField[r] += innerInterpolant[j] * externalField[edgeIndex[j]][r];
                }
            }
        }
        elemData.currentT = mesh.t;
        elemData.integralMonitor = mesh.integralMonitor;
        elemData.elemIndex = index;
        if (par.solutionAverage) {
            elemData.Wavg = calculateAvgW();
        }

        return u;
    }

    /**
     *
     * @param dt
     * @return L1norm(W - Wo)
     */
    public double calculateResiduumW(double dt
    ) {
        double rez = 0;
        for (int m = 0; m < nEqs; m++) {
            for (int j = 0; j < nBasis; j++) {
                rez = rez + Math.abs(W[m * nBasis + j] - Wo[m * nBasis + j]) / dt;
            }
        }

        return rez;
    }

    // ulozeni W do Wo
    public void copyW2Wo() {
        System.arraycopy(Wo, 0, Wo2, 0, nBasis * nEqs);
        System.arraycopy(W, 0, Wo, 0, nBasis * nEqs);
    }

    // ulozeni Wo do W (pri potizich s resenim)
    public void copyWo2W() {
        System.arraycopy(Wo, 0, W, 0, nBasis * nEqs);
    }

    //__________________________________________________________________________
    public double delta_t(double CFL) { //vypocet maximalniho casoveho kroku
        double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(centreVolumeInterpolant);
        double lam = eqn.maxEigenvalue(calculateAvgW(), elemData);
        double dt = CFL * elemSize / lam;

        return dt;
    }

    // vypocet geometrie _______________________________________________________
    public void computeGeometry() {
        Xes = new double[nFaces][];
        for (int k = 0; k < nFaces; k++) {
            Xes[k] = Int.faces[k].faceTransform.getXs();
        }
        S = new double[nFaces];
        for (int k = 0; k < nFaces; k++) {
            for (int i = 0; i < Int.faces[k].nIntEdge; i++) {
                S[k] = S[k] + Int.faces[k].JacobianFace[i] * Int.faces[k].weightsFace[i];
            }
        }

        n = new double[nFaces][][];
        for (int k = 0; k < nFaces; k++) {
            n[k] = new double[Int.faces[k].nIntEdge][dim];
            double[][] faceCoords = Int.faces[k].quadFace.getCoords();
            for (int p = 0; p < Int.faces[k].nIntEdge; p++) {
                n[k][p] = Int.faces[k].faceTransform.getNormal(faceCoords[p]);
            }
            //Mat.print(n[k][0]);
        }

        if (dim == 1) {  // only for 1D case
            n[0][0][0] = -1;
        }

        Xs = transform.getXs();

        area = 0;
        for (int i = 0; i < Int.nIntVolume; i++) {
            area = area + Int.JacobianVolume[i] * Int.weightsVolume[i];
        }
        //System.out.println(area);

        elemSize = 0;
        for (int k = 0; k < nFaces; k++) {
            elemSize += S[k];
        }
        elemSize = nFaces * area / elemSize;

        faceIndexReverse = new int[nFaces];
        for (int k = 0; k < nFaces; k++) {
            if (TT[k] > -1) {
                for (int s = 0; s < elems[TT[k]].nFaces; s++) {
                    if (elems[TT[k]].TT[s] == index) {
                        faceIndexReverse[k] = s;
                        break;
                    }
                }
            }
        }
    }

    public void nextTimeLevelMassMatrixes() {
        for (int i = 0; i < nBasis; i++) {
            for (int j = 0; j < nBasis; j++) {
                Mo2[i][j] = Mo[i][j];
                Mo[i][j] = M[i][j];
            }
        }
    }

    //____________________________________________________________________________
    public double[] calculateAvgW() { // funkce vraci stredni hodnotu W na kontrolnim objemu K
        double[] Ws = new double[nEqs];
        for (int j = 0; j < nEqs; j++) {
            for (int m = 0; m < nBasis; m++) {
                Ws[j] = Ws[j] + W[j * nBasis + m] * Is[m];
            }
        }
        return Ws;
    }

    public Equation getEqs() {
        return eqn;
    }

    public int getNEqs() {
        return nEqs;
    }

    public boolean geometryCheck(boolean writeOut) {
        boolean isOK = true;
        double[] test1 = new double[dim];
        double[] test2 = new double[dim];
        double sumTest1 = 0;
        double sumTest2 = 0;

        for (int k = 0; k < nFaces; k++) {
            Integration.Face face = Int.faces[k];
            for (int p = 0; p < face.nIntEdge; p++) { // edge integral
                double[] X = face.faceTransform.getX(face.quadFace.coords[p]);
                for (int d = 0; d < dim; d++) {
                    test1[d] += face.weightsFace[p] * face.JacobianFace[p] * n[k][p][d];
                    test2[d] += X[d] * face.weightsFace[p] * face.JacobianFace[p] * n[k][p][d];
                }
            }
        }
        for (int d = 0; d < dim; d++) {
            sumTest1 += test1[d];
            sumTest2 += test2[d] / dim;
            if (test2[d] < 0) {
                if (writeOut) {
                    System.out.println("Mesh element error, bad points order in file elements.txt!");
                }
                isOK = false;
            }
        }
        double testTol2 = 1e-11;
        if (elemType.order == 1) {
            testTol2 = 1e-4;
        }
        if (Math.abs(sumTest1) > 1e-11 || Math.abs(sumTest2 - area) > testTol2 || sumTest2 < 0) {
            if (writeOut) {
                System.out.println("Mesh element error! Control sum 1: " + Math.abs(sumTest1) + ", control sum 2: " + Math.abs(sumTest2 - area));
            }
            isOK = false;
        }

//            // check normals
//            boolean err = false;
//            for (int k = 0; k < nFaces; k++) {
//                Face face = Int.faces[k];
//                for (int p = 0; p < face.nIntEdge; p++) { // edge integral
//                    double[] X = face.faceTransform.getX(face.quadFace.coords[p]);
//                    double norOrient = Mat.scalar(n[k][p], Mat.minusVec(X, Xs));
//                    if (norOrient < 0) {
//                        err = true;
//
//                    }
//                }
//                if (err) {
//                    System.out.println("Mesh orientation error!");
//                }
//            }
        return isOK;
    }

    // optimalisation
    // for optimization toolbox, generate residuum(W)
    public void exportLocalFunctionalDerivative() {
        optimFunDer = new double[nBasis * nEqs];
        double[] V = new double[nBasis * nEqs];
        double Iw = optimalisationFunctional.combineFunctionals(computeFunctional(W));
        for (int i = 0; i < nBasis * nEqs; i++) {
            V[i] = par.h;
            double Iwh = optimalisationFunctional.combineFunctionals(computeFunctional(Mat.plusVec(W, V)));
            optimFunDer[i] = (Iwh - Iw) / par.h;
            V[i] = 0;
        }
    }

    // for optimization toolbox, generate residuum(W)
    public void exportLocalR() {
        if (ti.isImplicit) {
            double[] V = new double[nBasis * nEqs];
            double[] Rw = new double[nBasis * nEqs];
            double[][] RwN = new double[nFaces][];
            residuum(V, Rw, RwN);
            System.arraycopy(Rw, 0, ((Implicit) ti).RHS_loc, 0, Rw.length);
        } else {
            throw new UnsupportedOperationException("operation not supported for this time integration method");
        }
    }

    // for optimization toolbox, generate only residuum derivation (dR/dW)
    public void exportLocalJacobiMatrix() {
        if (ti.isImplicit) {
            // vnitrni element - krivkovy i objemovy integral
            double[] V = new double[nBasis * nEqs];
            double[] Rw = new double[nBasis * nEqs];
            double[][] RwNeigh = new double[nFaces][];
            double[][] RwNeighH = new double[nFaces][];
            for (int k = 0; k < nFaces; k++) {
                if (TT[k] > -1) {
                    RwNeigh[k] = new double[elems[TT[k]].nBasis * nEqs];
                    RwNeighH[k] = new double[elems[TT[k]].nBasis * nEqs];
                }
            }
            // compute residuum
            residuum(V, Rw, RwNeigh);
            double[][] ADiag = ((Implicit) ti).ADiag;
            double h = par.h;
            for (int i = 0; i < nBasis * nEqs; i++) {
                for (int j = 0; j < Rw.length; j++) {
                    ADiag[i][j] = -Rw[j];
                }
                V[i] = h;
                residuum(V, ADiag[i], RwNeighH);
                V[i] = 0;
                for (int k = 0; k < nFaces; k++) {
                    if (TT[k] > -1) {
                        double[][] Aaux = ((Implicit) elems[TT[k]].ti).ANeighs[faceIndexReverse[k]].A;
                        for (int j = 0; j < RwNeighH[k].length; j++) {
                            Aaux[i][j] = (RwNeighH[k][j] - RwNeigh[k][j]) / h;
                        }
                    }
                }
            }
            Mat.divide(ADiag, -h);
        } else {
            throw new UnsupportedOperationException("operation not supported for this time integration method");
        }
    }

    public double[] computeFunctional(double[] V) {
        if (V == null) {
            V = new double[nBasis * nEqs];
        }
        int nF = optimalisationFunctional.getN();
        double[] f = new double[nF];

        // vypocet toku hranici
        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
            Mat.plusVecToVec(f, computeFunctionalsWall(k, V));
        }

        for (int p = 0; p < Int.nIntVolume; p++) {
            double[] base = Int.basisVolume[p];
            double[][] dBase = Int.dXbasisVolume[p];
            double Jac = Int.JacobianVolume[p];
            double weight = Int.weightsVolume[p];

            // interpolation of mesh velocity
            double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[p]);

            double[] WInt = new double[nEqs];
            double[] dWInt = new double[dim * nEqs];
            for (int m = 0; m < nEqs; m++) {
                for (int j = 0; j < nBasis; j++) {
                    WInt[m] += (W[m * nBasis + j] + V[m * nBasis + j]) * base[j];
                    for (int d = 0; d < dim; d++) {
                        dWInt[nEqs * d + m] += (W[m * nBasis + j] + V[m * nBasis + j]) * dBase[j][d];
                    }
                }
            }
            double[] aux = optimalisationFunctional.insideValue(WInt, dWInt, elemData);
            for (int i = 0; i < nF; i++) {
                f[i] += Jac * weight * aux[i];
            }
        }
        return f;
    }

    // tato funkce vypocitava reziduum__________________________________________
    public double[] computeFunctionalsWall(int k, double[] V) {
        int nF = optimalisationFunctional.getN();
        double[] f = new double[nF];
        int[] edgeIndex = Int.faces[k].faceIndexes;

        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
            double Jac = Int.faces[k].JacobianFace[p];
            double weight = Int.faces[k].weightsFace[p];

            // interpolation of mesh velocity
            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);

            double[] WL = new double[nEqs];
            double[] dWL = new double[dim * nEqs];

            // values from boundary inlet (WL, dWL)
            for (int j = 0; j < nEqs; j++) {
                for (int m = 0; m < nBasis; m++) {
                    WL[j] += (W[j * nBasis + m] + V[j * nBasis + m]) * baseLeft[m];
                    for (int d = 0; d < dim; d++) {
                        dWL[nEqs * d + j] += (W[j * nBasis + m] + V[j * nBasis + m]) * dBaseLeft[m][d];
                    }
                }
            }
            double[] aux = optimalisationFunctional.boundaryValue(WL, dWL, n[k][p], TT[k], elemData);
            for (int i = 0; i < nF; i++) {
                f[i] += Jac * weight * aux[i];
            }
        }
        return f;
    }

    public double[] computeIntegralSolutionMonitor(int nIntegralMonitor) {
        double[] f = new double[nIntegralMonitor];

        // vypocet toku hranici
        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
            computeIntegralSolutionMonitorOnWall(k, f);
        }

        for (int p = 0; p < Int.nIntVolume; p++) {
            double[] base = Int.basisVolume[p];
            double[][] dBase = Int.dXbasisVolume[p];
            double Jac = Int.JacobianVolume[p];
            double weight = Int.weightsVolume[p];

            // interpolation of mesh velocity
            double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[p]);

            double[] WInt = new double[nEqs];
            double[] dWInt = new double[dim * nEqs];
            for (int m = 0; m < nEqs; m++) {
                for (int j = 0; j < nBasis; j++) {
                    WInt[m] += W[m * nBasis + j] * base[j];
                    for (int d = 0; d < dim; d++) {
                        dWInt[nEqs * d + m] += W[m * nBasis + j] * dBase[j][d];
                    }
                }
            }
            double[] aux = mesh.solMonitor.insideValue(WInt, dWInt, elemData);
            if (aux != null) {
                for (int r = 0; r < aux.length; r++) {
                    f[r] += Jac * weight * aux[r];
                }
            }
        }
        return f;
    }

    // tato funkce vypocitava reziduum__________________________________________
    public void computeIntegralSolutionMonitorOnWall(int k, double[] f) {
        int[] edgeIndex = Int.faces[k].faceIndexes;

        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
            double Jac = Int.faces[k].JacobianFace[p];
            double weight = Int.faces[k].weightsFace[p];

            // interpolation of mesh velocity
            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);

            double[] WL = new double[nEqs];
            double[] dWL = new double[dim * nEqs];

            // values from boundary inlet (WL, dWL)
            for (int j = 0; j < nEqs; j++) {
                for (int m = 0; m < nBasis; m++) {
                    WL[j] += W[j * nBasis + m] * baseLeft[m];
                    for (int d = 0; d < dim; d++) {
                        dWL[nEqs * d + j] += W[j * nBasis + m] * dBaseLeft[m][d];
                    }
                }
            }
            double[] aux = mesh.solMonitor.boundaryValue(WL, dWL, n[k][p], TT[k], elemData);
            if (aux != null) {
                for (int r = 0; r < aux.length; r++) {
                    f[r] += Jac * weight * aux[r];
                }
            }
        }
    }

    public TimeIntegrationElement getTimeIntegrationElement(String methodClassName, Element elem) throws IOException {

        String className = "flowpro.core.element." + methodClassName;

        try {
            Class<TimeIntegrationElement> tiClass = (Class<TimeIntegrationElement>) Class.forName("flowpro.core.element." + methodClassName);
            TimeIntegrationElement tiElem = (TimeIntegrationElement) tiClass.newInstance();
            tiElem.set(elem);
            return tiElem;
        } catch (ClassNotFoundException ex) {
            throw new IOException("class \"" + className + "\" not found", ex);
        } catch (InstantiationException | IllegalAccessException | SecurityException | IllegalArgumentException ex) {
            throw new IOException("error while loading class \"" + className + "\": " + ex, ex);
        }
    }
}
