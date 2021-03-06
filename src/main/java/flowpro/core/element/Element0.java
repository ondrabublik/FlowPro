///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//
//package flowpro.core.element;
//
//import flowpro.api.ElementData;
//import flowpro.api.Equation;
//import flowpro.api.Functional;
//import flowpro.api.Mat;
//import flowpro.core.Integration;
//import flowpro.core.Mesh;
//import flowpro.core.Parameters;
//import flowpro.core.basis.Basis;
//import flowpro.core.curvedBoundary.FaceCurvature;
//import flowpro.core.elementType.ElementType;
//import flowpro.core.quadrature.QuadratureCentral;
//import flowpro.core.transformation.Transformation;
//import java.io.IOException;
//import java.io.Serializable;
//import java.util.Arrays;
//
///**
// *
// * @author obublik
// */
///**
// * Represents an element in the mesh.
// */
//
//public class Element0 implements Serializable {
//
//    // mesh
//    private final Mesh mesh;
//    // parameters
//    Parameters par;
//    // equations
//    int nEqs;
//    Equation eqn;
//    // mesh elements
//    Element[] elems;
//    
//    /* constants */
//    public double c_IP;  // ???
//    public double eps;  // ???
//
//    /* geometry */
//    public final int[] TT;     // indexes of neighbours
//    public final int[] TP;
//    public final int[] TEale;     // type of boundary for ALE
//    public final int[] TEshift;
//    public double[][] shift;
//    public FaceCurvature fCurv;
//    public int nFaces;       // number of faces
//    public final int nVertices;       // number of vertices
//    public double[][] vertices;  // components of vertices [xi,yi,zi]
//    public double[][] verticesOld;
//    public double[][] verticesOld2;
//    public double[][] vertices0;  // components of vertices [xi,yi,zi]
//    public double[][] U;  // x-components of velocity [ui, vi, wi]
//    public double[][] Uinit;
//    public double[][][] n;  // x-components of normals
//    public double[] S;
//    public double area;  // element area
//    public double[] wallDistance; // vzdalennost od steny
//    public double[] Xs;  // coordinates of element's midpoint
//    public double[][] Xes; // coordinates of element's edge midpoint
//    // size of element defined as 2 times the diameter of inscribed cyrcle
//    public double elemSize;
//    public boolean insideComputeDomain; // true if element is insideComputeDomain subdomain, false on the boundary
//    public boolean insideMetisDomain; // true insideComputeDomain computational subdomain (producet by metis)
//    public boolean insideAssemblerDomain; // interior for Jacobi matrix assemble
//    public boolean gmresLoad, gmresSave;
//    public double[][] blendFun; // blending function used for calculation of new mesh position
//
//    // mesh
//    public final int index;
//    public int[] faceIndexReverse;
//
//    public final int dim;
//    public int nBasis; 	  // pocet bazovych funkci
//    public Transformation transform;
//    public Basis basis;
//    public Integration Int;
//    QuadratureCentral qRules;
//    public ElementType elemType;
//
//    //private double k_el; // kvalita elementu
//    public double[] Is; // integral bazove funkce
//    public double[][] M;	 // matice hmotnosti
//    double[][] iM; // inverze matice hmotnosti (pouze pro implicitni metodu)
//    double[][] Mo;	 // matice hmotnosti v predchozi casove hladine
//    double[][] Mo2;	 // matice hmotnosti
//    public double[][] ADiag;
//    double[][] PrecondJacobi; // inverze ADiag
//    public double[] RHS_loc;
//    public Neighbour[] ANeighs; //matice sousedu
//
//    // damping
//    double[][] TrunOrd;
//    public double[] dampInner;
//    public double[] innerIndicatorOld;
//
//    // local time step
//    public double tLTS;
//    public double tLTSold;
//    public int stepLTS;
//    public double[] WLTS, W1LTS, W1LTSo, W2LTS, W2LTSo;   // auxilary for local time stepping
//
//    // parametery pro prenos do equations
//    public double[] centreVolumeInterpolant;
//    public final ElementData elemData;
//
//    // optimalisation
//    public Functional optimalisationFunctional; // functional for optimalisation
//    double[] optimFunDer;
//
//    // matice globalnich indexu a globalni matice s pravou stranou
//    public int[] gi_U;
//    
//    double[] initW;
//    public double[] W;     // hodnota W v n+1 te casove hladine
//    public double[] Wo;    // hodnota v n te casove hladine
//    public double[] Wo2;   // hodnota v n-1 casove hladine
//
//    // derivative correction
//    double[] R;
//
//    // external field
//    double[][] externalField;
//
//    public Element(int index, double[][] vertices, double[][] Uinit, double[] wallDistance, double[][] externalField, int[] TT, int[] TP, int[] TEale, int[] TEshift, double[][] shift, FaceCurvature fCurv, double[][] blendFun, double[] initW,
//            Mesh mesh, ElementType elemType) throws IOException {
//        
//        this.mesh = mesh;
//        eqn = mesh.getEqn();
//        dim = eqn.dim();
//        nEqs = eqn.nEqs();
//        par = mesh.getPar();
//        qRules = mesh.getQRules();
//        elems = mesh.getElems();
//        
//        this.index = index;
//        this.TT = TT;
//        this.TP = TP;
//        this.TEale = TEale;
//        this.TEshift = TEshift;
//        this.shift = shift;
//        this.fCurv = fCurv;
//        this.vertices = vertices;
//        this.Uinit = Uinit;
//        this.blendFun = blendFun;
//        this.wallDistance = wallDistance;
//        this.externalField = externalField;
//        this.nFaces = elemType.nFaces;
//        this.initW = initW;
//        this.nVertices = vertices.length;
//        this.elemType = elemType;
//        
//
//        insideComputeDomain = true;
//        gmresLoad = false;
//        gmresSave = false;
//
//        U = new double[nVertices][dim];
//        verticesOld = new double[nVertices][dim];
//        verticesOld2 = new double[nVertices][dim];
//
//        elemData = new ElementData(dim);
//        centreVolumeInterpolant = new double[nVertices];
//        Arrays.fill(centreVolumeInterpolant, 1.0 / nVertices);
//
//        // damping
//        dampInner = new double[nEqs];
//        innerIndicatorOld = new double[nEqs];
//    }
//
//    public void initBasis() throws IOException {
//        transform = elemType.getVolumeTransformation(vertices, fCurv, par);
//        basis = elemType.getBasis(transform);
//        nBasis = basis.nBasis;
//    }
//
//    public void initIntegration() throws IOException {
//        Int = new Integration(elemType, dim, basis, transform, TT, TEshift, shift, qRules);
//
//        if (!par.isExplicit) {
//            ADiag = new double[nEqs * nBasis][nEqs * nBasis];
//            RHS_loc = new double[nEqs * nBasis];
//        }
//    }
//
//    public void initCondition() throws IOException {
//        vertices0 = new double[nVertices][dim];
//        for (int d = 0; d < dim; d++) {
//            for (int i = 0; i < nVertices; i++) {
//                vertices0[i][d] = vertices[i][d];
//                U[i][d] = Uinit[i][d];
//            }
//        }
//
//        // fill the solution vector with initial condition
//        W = new double[nBasis * nEqs];
//        Wo = new double[nBasis * nEqs];
//        Wo2 = new double[nBasis * nEqs];
//        if (initW.length == nBasis * nEqs) {
//            System.arraycopy(initW, 0, W, 0, nBasis * nEqs);
//            System.arraycopy(initW, 0, Wo, 0, nBasis * nEqs);
//            System.arraycopy(initW, 0, Wo2, 0, nBasis * nEqs);
//        } else {
//            switch (basis.basisType) {
//                case "lagrange":
//                    for (int i = 0; i < nBasis; i++) {
//                        for (int j = 0; j < nEqs; j++) {
//                            W[j * nBasis + i] = initW[j];
//                            Wo[j * nBasis + i] = initW[j];
//                            Wo2[j * nBasis + i] = initW[j];
//                        }
//                    }
//                    break;
//                case "orthogonal":
//                case "taylor":
//                    for (int j = 0; j < nEqs; j++) {
//                        W[j * nBasis] = initW[j];
//                        Wo[j * nBasis] = initW[j];
//                        Wo2[j * nBasis] = initW[j];
//                    }
//                    break;
//            }
//        }
//    }
//
//    public void computeInvertMassMatrix() {
//        // vypocitava inverzi diagonaly, ktera se pouzije pro predpodminovac
//        double[][] invM = Mat.invert(M);
//        iM = new double[nEqs * nBasis][nEqs * nBasis];
//        for (int m = 0; m < nEqs; m++) {
//            for (int i = 0; i < nBasis; i++) {
//                for (int j = 0; j < nBasis; j++) {
//                    iM[m * nBasis + i][m * nBasis + j] = invM[i][j];
//                }
//            }
//        }
//    }
//
//    // Aplikace limiteru
//    public void limiter(boolean isFirstIter) {
//        eps = 0;
//        c_IP = 0;
//        if (elemType.order > 1) {
//            // max eigenvalue
//            double lam = eqn.maxEigenvalue(calculateAvgW(), elemData);
//
//            // artificial viscosity
//            if (par.dampTol > 0) {
//                double g_shock = shock_senzor(par.dampTol);
//                eps = lam * elemSize / elemType.order * g_shock;
//            }
//
//            if (par.dampInnerTol > 0) {
//                // inner damping based on artificial viscosity
//                double[] g_shockInner = eqn.combineShockSensors(shock_senzor_all_eqn(par.dampInnerTol));
//                if (isFirstIter && !par.continueComputation) {
//                    for (int i = 0; i < g_shockInner.length; i++) {
//                        // g_shockInner[i] = 1;
//                    }
//                }
//
//                for (int m = 0; m < nEqs; m++) {
//                    double indicator = Math.max(g_shockInner[m], innerIndicatorOld[m]);
//                    innerIndicatorOld[m] = 0.1 * indicator;
//                    if (par.dampInnerCoef != null) {
//                        if (par.dampInnerCoef.length == 1) {
//                            dampInner[m] = par.dampInnerCoef[0] * indicator;
//                        } else {
//                            dampInner[m] = par.dampInnerCoef[m] * indicator;
//                        }
//                    } else {
//                        dampInner[m] = lam * elemSize / elemType.order * indicator;
//                    }
//                }
//            }
//
//            if (eqn.isDiffusive()) {
//                c_IP = par.penalty; // * elemSize;
//            }
//        }
//    }
//
//    public void computeExplicitStep(double dtStep) {
//        double[] V = new double[nBasis * nEqs];
//        double[] Rw = new double[nBasis * nEqs];
//        double[][] RwN = new double[nFaces][];
//        limiter(false);
//        switch (par.orderInTime) {
//            case 2: // RK2
//                System.arraycopy(W, 0, WLTS, 0, nBasis * nEqs);
//                stepLTS = 0;
//                residuum(V, Rw, RwN);
//                Rw = Mat.times(iM, Rw);
//                for (int j = 0; j < W.length; j++) {
//                    W1LTSo[j] = W1LTS[j];
//                    W[j] = WLTS[j] + dtStep / 2 * Rw[j];
//                    W1LTS[j] = W[j];
//                }
//
//                stepLTS = 1;
//                Arrays.fill(Rw, 0);
//                residuum(V, Rw, RwN);
//                Rw = Mat.times(iM, Rw);
//                for (int j = 0; j < W.length; j++) {
//                    W[j] = WLTS[j] + dtStep * Rw[j];
//                }
//                break;
//            case 3: // RK3
//                System.arraycopy(W, 0, WLTS, 0, nBasis * nEqs);
//                stepLTS = 0;
//                residuum(V, Rw, RwN);
//                Rw = Mat.times(iM, Rw);
//                for (int j = 0; j < W.length; j++) {
//                    W1LTSo[j] = W1LTS[j];
//                    W[j] = WLTS[j] + dtStep * Rw[j];
//                    W1LTS[j] = W[j];
//                }
//
//                stepLTS = 1;
//                Arrays.fill(Rw, 0);
//                residuum(V, Rw, RwN);
//                Rw = Mat.times(iM, Rw);
//                for (int j = 0; j < W.length; j++) {
//                    W2LTSo[j] = W2LTS[j];
//                    W[j] = 0.75 * WLTS[j] + 0.25 * (W1LTS[j] + dtStep * Rw[j]);
//                    W2LTS[j] = W[j];
//                }
//
//                stepLTS = 2;
//                Arrays.fill(Rw, 0);
//                residuum(V, Rw, RwN);
//                Rw = Mat.times(iM, Rw);
//                for (int j = 0; j < W.length; j++) {
//                    W[j] = 1.0 / 3 * WLTS[j] + 2.0 / 3 * (W2LTS[j] + dtStep * Rw[j]);
//                }
//        }
//        //eqn.limitUnphysicalValues(calculateAvgW(), W, nBasis);
//    }
//
//    void assembleRHS(double[] Rw, double[] a1, double[] a2, double[] a3) {
//        System.arraycopy(Rw, 0, RHS_loc, 0, Rw.length);
//        for (int m = 0; m < nEqs; m++) {
//            for (int i = 0; i < nBasis; i++) {
//                for (int j = 0; j < nBasis; j++) {
//                    RHS_loc[nBasis * m + i] -= M[i][j] * a1[m] * W[m * nBasis + j] + Mo[i][j] * a2[m] * Wo[m * nBasis + j] + Mo2[i][j] * a3[m] * Wo2[m * nBasis + j];
//                }
//            }
//        }
//    }
//
//    // Generovani radku globalni matice a vektoru prave strany
//    public void assembleJacobiMatrix(double[] a1, double[] a2, double[] a3, double[] dual) {
//
//        // vnitrni element - krivkovy i objemovy integral
//        double[] V = new double[nBasis * nEqs];
//        double[] Rw = new double[nBasis * nEqs];
//        double[][] RwNeigh = new double[nFaces][];
//        double[][] RwNeighH = new double[nFaces][];
//        for (int k = 0; k < nFaces; k++) {
//            if (TT[k] > -1) {
//                RwNeigh[k] = new double[elems[TT[k]].nBasis * nEqs];
//                RwNeighH[k] = new double[elems[TT[k]].nBasis * nEqs];
//            }
//        }
//        // compute residuum
//        residuum(V, Rw, RwNeigh);
//        // assemble rhs
//        assembleRHS(Rw, a1, a2, a3);
//
//        if (par.useJacobiMatrix && eqn.isEquationsJacobian()) { // fast assemble when jacobian of equations is known
//            residuumWithJacobian(ADiag, ANeighs);
//        } else { // slow assemble when jacobian of equations is unknown
//            double h = par.h;
//            for (int i = 0; i < nBasis * nEqs; i++) {
//                for (int j = 0; j < Rw.length; j++) {
//                    ADiag[i][j] = -Rw[j];
//                }
//                V[i] = h;
//                residuum(V, ADiag[i], RwNeighH);
//                V[i] = 0;
//                for (int k = 0; k < nFaces; k++) {
//                    if (TT[k] > -1 && elems[TT[k]].insideComputeDomain) {
//                        double[][] Aaux = elems[TT[k]].ANeighs[faceIndexReverse[k]].A;
//                        for (int j = 0; j < RwNeighH[k].length; j++) {
//                            Aaux[i][j] = (RwNeighH[k][j] - RwNeigh[k][j]) / h;
//                        }
//                    }
//                }
//            }
//            Mat.divide(ADiag, -h);
//        }
//
//        // pricteni matice hmotnosti
//        for (int m = 0; m < nEqs; m++) {
//            for (int i = 0; i < nBasis; i++) {
//                for (int j = 0; j < nBasis; j++) {
//                    ADiag[nBasis * m + i][nBasis * m + j] += (a1[m] + dual[m]) * M[i][j];
//                }
//            }
//        }
//    }
//
//    private void residuumWithJacobian(double[][] ADiag, Neighbour[] Sous) {
//        // vypocet toku hranici
//        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
//            residuumWithJacobianWall(k, ADiag, Sous[k]);
//        }
//
//        if (elemType.order > 1) { // volume integral only for DGFEM
//            double[] nor = new double[dim];
//            double[] a = null;
//            double[] ad = null;
//            double[] ap = null;
//
//            for (int p = 0; p < Int.nIntVolume; p++) {
//                double[] base = Int.basisVolume[p];
//                double[][] dBase = Int.dXbasisVolume[p];
//                double Jac = Int.JacobianVolume[p];
//                double weight = Int.weightsVolume[p];
//
//                // interpolation of mesh velocity, and other data
//                double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[p]);
//
//                double[] WInt = new double[nEqs];
//                double[] dWInt = new double[dim * nEqs];
//                for (int m = 0; m < nEqs; m++) {
//                    for (int j = 0; j < nBasis; j++) {
//                        WInt[m] += W[m * nBasis + j] * base[j];
//                        for (int d = 0; d < dim; d++) {
//                            dWInt[nEqs * d + m] += W[m * nBasis + j] * dBase[j][d];
//                        }
//                    }
//                }
//                // convection
//                for (int d = 0; d < dim; d++) {
//                    if (eqn.isConvective()) {
//                        nor[d] = 1;
//                        a = eqn.convectiveFluxJacobian(WInt, nor, elemData);
//                        nor[d] = 0;
//                        for (int m = 0; m < nEqs; m++) {
//                            a[nEqs * m + m] -= u[d];
//                        }
//                    }
//                    if (eqn.isDiffusive()) {
//                        nor[d] = 1;
//                        ad = eqn.diffusiveFluxJacobian(WInt, dWInt, nor, elemData);
//                        nor[d] = 0;
//                    }
//                    if (eqn.isSourcePresent()) {
//                        nor[d] = 1;
//                        ap = eqn.sourceTermJacobian(WInt, dWInt, elemData);
//                        nor[d] = 0;
//                    }
//                    for (int m = 0; m < nEqs; m++) {
//                        for (int q = 0; q < nEqs; q++) {
//                            for (int i = 0; i < nBasis; i++) {
//                                for (int j = 0; j < nBasis; j++) {
//                                    if (eqn.isConvective()) {
//                                        ADiag[nBasis * m + i][nBasis * q + j] -= Jac * weight * a[nEqs * q + m] * base[i] * dBase[j][d];
//                                        if (m == q) {
//                                            ADiag[nBasis * m + i][nBasis * q + j] += (eps + par.dampConst) * Jac * weight * dBase[i][d] * dBase[j][d];
//                                        }
//                                    }
//                                    if (eqn.isDiffusive()) {
//                                        ADiag[nBasis * m + i][nBasis * q + j] += Jac * weight * ad[nEqs * q + m] * base[i] * dBase[j][d];
//                                        for (int r = 0; r < dim; r++) {
//                                            ADiag[nBasis * m + i][nBasis * q + j] += Jac * weight * ad[nEqs * nEqs * (r + 1) + nEqs * q + m] * dBase[i][r] * dBase[j][d];
//                                        }
//                                    }
//                                    if (eqn.isSourcePresent()) {
//                                        ADiag[nBasis * m + i][nBasis * q + j] -= Jac * weight * ap[nEqs * q + m] * base[i] * dBase[j][d] + Jac * weight * ap[nEqs * nEqs * (d + 1) + nEqs * q + m] * dBase[i][d] * dBase[j][d];
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    private void residuumWithJacobianWall(int k, double[][] ADiag, Neighbour Sous) {
//        double[] aL = null;
//        double[] aR = null;
//        double[] adL = null;
//        double[] adR = null;
//        int[] edgeIndex = Int.faces[k].faceIndexes;
//        double h = 1e-6;
//        double[] V = new double[nEqs];
//
//        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
//            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
//            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
//            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
//            double Jac = Int.faces[k].JacobianFace[p];
//            double weight = Int.faces[k].weightsFace[p];
//            double[] baseRight = null;
//            double[][] dBaseRight = null;
//            if (TT[k] > -1) {
//                baseRight = Int.faces[k].basisFaceRight[p];
//                dBaseRight = Int.faces[k].dXbasisFaceRight[p];
//            }
//
//            // interpolation of mesh velocity
//            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);
//
//            double dL = 0;
//            for (int d = 0; d < dim; d++) {
//                if (TT[k] > -1) {
//                    dL += (Xs[d] - elems[TT[k]].Xs[d]) * n[k][p][d];
//                } else {
//                    dL += (Xs[d] - Xes[k][d]) * n[k][p][d];
//                }
//            }
//            dL = Math.abs(dL);
//
//            double[] WL = new double[nEqs];
//            double[] WR = new double[nEqs];
//            double[] dWL = new double[dim * nEqs];
//            double[] dWR = new double[dim * nEqs];
//
//            // values from boundary inlet (WL, dWL)
//            for (int j = 0; j < nEqs; j++) {
//                for (int m = 0; m < nBasis; m++) {
//                    WL[j] += W[j * nBasis + m] * baseLeft[m];
//                    for (int d = 0; d < dim; d++) {
//                        dWL[nEqs * d + j] += W[j * nBasis + m] * dBaseLeft[m][d];
//                    }
//                }
//            }
//
//            // values from boundary outlet (WR, dWR)
//            if (TT[k] > -1) {
//                double[] WRp = elems[TT[k]].getW(par.isExplicit, tLTS);
//                for (int m = 0; m < nEqs; m++) {
//                    int nRBasis = elems[TT[k]].nBasis;
//                    for (int j = 0; j < nRBasis; j++) {
//                        WR[m] += WRp[m * nRBasis + j] * baseRight[j];
//                        for (int d = 0; d < dim; d++) {
//                            dWR[nEqs * d + m] += WRp[m * nRBasis + j] * dBaseRight[j][d];
//                        }
//                    }
//                }
//            } else {
//                WR = eqn.boundaryValue(WL, n[k][p], TT[k], elemData);
//                System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
//            }
//
//            // inviscid flux in integration point
//            if (eqn.isConvective()) {
//                double vn = 0;
//                for (int d = 0; d < dim; d++) {
//                    vn = vn + u[d] * n[k][p][d];
//                }
//                aL = eqn.convectiveFluxJacobian(WL, n[k][p], elemData);	// nevazky tok
//                if (TT[k] > -1) {
//                    aR = eqn.convectiveFluxJacobian(WR, n[k][p], elemData);	// nevazky tok
//                }
//            }
//
//            // viscouse flux in integration point
//            if (eqn.isDiffusive()) {
//                adL = eqn.diffusiveFluxJacobian(WL, dWL, n[k][p], elemData);	// vazky tok
//                if (TT[k] > -1) {
//                    adR = eqn.diffusiveFluxJacobian(WR, dWR, n[k][p], elemData);	// vazky tok
//                }
//            }
//
//            if (TT[k] > -1) { // inner edge
//                int nRBasis = elems[TT[k]].nBasis;
//                double[] dBazeSumL = new double[nBasis];
//                if (TT[k] > -1) {
//                    for (int i = 0; i < nBasis; i++) {
//                        for (int d = 0; d < dim; d++) {
//                            dBazeSumL[i] += dBaseLeft[i][d] * n[k][p][d];
//                        }
//                    }
//                }
//                double[] dBazeSumR = new double[nRBasis];
//                if (TT[k] > -1) {
//                    for (int i = 0; i < nRBasis; i++) {
//                        for (int d = 0; d < dim; d++) {
//                            dBazeSumR[i] += dBaseRight[i][d] * n[k][p][d];
//                        }
//                    }
//                }
//
//                // LF schema ==================================
//                if (eqn.isConvective()) {
//                    double[] Ws = new double[nEqs];
//                    for (int m = 0; m < nEqs; m++) {
//                        Ws[m] = (WL[m] + WR[m]) / 2;
//                    }
//                    double lam = eqn.maxEigenvalue(Ws, elemData);
//                    for (int m = 0; m < nEqs; m++) {
//                        aL[nEqs * m + m] += lam;
//                        aR[nEqs * m + m] -= lam;
//                    }
//
//                    for (int m = 0; m < nEqs; m++) {
//                        V[m] = h;
//                        double lamh = eqn.maxEigenvalue(Mat.plusVec(Ws, V), elemData);
//                        double dlam = (lamh - lam) / h;
//                        V[m] = 0;
//                        for (int q = 0; q < nEqs; q++) {
//                            aL[nEqs * q + m] += dlam * WL[q];
//                            aR[nEqs * q + m] -= dlam * WR[q];
//                        }
//                    }
//                }
//                // end LF schema ==================================
//
//                // IP =============================================
//                if (eqn.isDiffusive()) {
//                    for (int m = 0; m < nEqs; m++) {
//                        adL[nEqs * m + m] -= (c_IP + elems[TT[k]].c_IP) / 2;
//                        adR[nEqs * m + m] += (c_IP + elems[TT[k]].c_IP) / 2;
//                    }
//                }
//                // end IP =========================================
//
//                for (int m = 0; m < nEqs; m++) {
//                    for (int q = 0; q < nEqs; q++) {
//                        for (int i = 0; i < nBasis; i++) {
//                            for (int j = 0; j < nBasis; j++) {
//                                if (eqn.isConvective()) {
//                                    ADiag[nBasis * m + i][nBasis * q + j] += 0.5 * Jac * weight * aL[nEqs * q + m] * baseLeft[i] * baseLeft[j];
//                                    if (m == q) {
//                                        ADiag[nBasis * m + i][nBasis * q + j] -= 0.5 * (0.5 * (eps + elems[TT[k]].eps) + par.dampConst) * Jac * weight * dBazeSumL[i] * baseLeft[j];
//                                    }
//                                }
//                                if (eqn.isDiffusive()) {
//                                    ADiag[nBasis * m + i][nBasis * q + j] -= 0.5 * Jac * weight * adL[nEqs * q + m] * baseLeft[i] * baseLeft[j];
//                                    for (int d = 0; d < dim; d++) {
//                                        ADiag[nBasis * m + i][nBasis * q + j] -= 0.5 * Jac * weight * adL[nEqs * nEqs * (d + 1) + nEqs * q + m] * dBaseLeft[i][d] * baseLeft[j];
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//
//                for (int m = 0; m < nEqs; m++) {
//                    for (int q = 0; q < nEqs; q++) {
//                        for (int i = 0; i < nRBasis; i++) {
//                            for (int j = 0; j < nBasis; j++) {
//                                if (eqn.isConvective()) {
//                                    Sous.A[nRBasis * m + i][nBasis * q + j] += 0.5 * Jac * weight * aR[nEqs * q + m] * baseRight[i] * baseLeft[j];
//                                    if (m == q) {
//                                        Sous.A[nRBasis * m + i][nBasis * q + j] -= 0.5 * (0.5 * (eps + elems[TT[k]].eps) + par.dampConst) * Jac * weight * dBazeSumR[i] * baseLeft[j];
//                                    }
//                                }
//                                if (eqn.isDiffusive()) {
//                                    Sous.A[nRBasis * m + i][nBasis * q + j] -= 0.5 * Jac * weight * adR[nEqs * q + m] * baseRight[i] * baseLeft[j];
//                                    for (int d = 0; d < dim; d++) {
//                                        Sous.A[nRBasis * m + i][nBasis * q + j] -= 0.5 * Jac * weight * adR[nEqs * nEqs * (d + 1) + nEqs * q + m] * dBaseRight[i][d] * baseLeft[j];
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            } else { // boundary edge
//                if (eqn.isConvective()) {
//                    double[] fR0 = eqn.numericalConvectiveFlux(WL, WR, n[k][p], TT[k], elemData);
//                    for (int m = 0; m < nEqs; m++) {
//                        V[m] = h;
//                        double[] WLh = Mat.plusVec(WL, V);
//                        double[] WRh = eqn.boundaryValue(WLh, n[k][p], TT[k], elemData);
//                        double[] fRh = eqn.numericalConvectiveFlux(WLh, WRh, n[k][p], TT[k], elemData);
//                        V[m] = 0;
//                        for (int q = 0; q < nEqs; q++) {
//                            double derfR = (fRh[q] - fR0[q]) / h;
//                            for (int i = 0; i < nBasis; i++) {
//                                for (int j = 0; j < nBasis; j++) {
//                                    ADiag[nBasis * m + i][nBasis * q + j] += Jac * weight * derfR * baseLeft[i] * baseLeft[j];
//                                }
//                            }
//                        }
//                    }
//                } else {
//                    for (int m = 0; m < nEqs; m++) {
//                        for (int q = 0; q < nEqs; q++) {
//                            for (int i = 0; i < nBasis; i++) {
//                                for (int j = 0; j < nBasis; j++) {
//                                    ADiag[nBasis * m + i][nBasis * q + j] += Jac * weight * aL[nEqs * q + m] * baseLeft[i] * baseLeft[j];
//                                }
//                            }
//                        }
//                    }
//                }
//
//                if (eqn.isDiffusive()) {
//                    double beta0 = par.beta0;
//                    if (par.order == 1) {
//                        beta0 = 0;
//                    }
//                    double[] Wc = new double[nEqs];
//                    double[] dWc = new double[nEqs * dim];
//                    for (int m = 0; m < nEqs; m++) {
//                        if (TT[k] > 0) {
//                            Wc[m] = (WL[m] + WR[m]) / 2;
//                        } else {
//                            Wc[m] = WR[m];
//                        }
//                        for (int d = 0; d < dim; d++) {
//                            dWc[nEqs * d + m] = (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2 + beta0 * (WR[m] - WL[m]) / dL * n[k][p][d];
//                        }
//                    }
//                    double[] fvR0 = eqn.numericalDiffusiveFlux(Wc, dWc, n[k][p], TT[k], elemData);
//                    for (int m = 0; m < nEqs; m++) {
//                        V[m] = h;
//                        double[] WLh = Mat.plusVec(WL, V);
//                        double[] WRh = eqn.boundaryValue(WLh, n[k][p], TT[k], elemData);
//                        double[] Wch = new double[nEqs];
//                        for (int q = 0; q < nEqs; q++) {
//                            if (TT[k] > 0) {
//                                Wch[q] = (WLh[q] + WRh[q]) / 2;
//                            } else {
//                                Wch[q] = WRh[q];
//                            }
//                        }
//                        double[] fvRh = eqn.numericalDiffusiveFlux(Wch, dWc, n[k][p], TT[k], elemData);
//                        for (int q = 0; q < nEqs; q++) {
//                            double derfvR = (fvRh[q] - fvR0[q]) / h;
//                            double derIP = c_IP * (WLh[q] - WL[q] - (WRh[q] - WR[q])) / h;
//                            for (int i = 0; i < nBasis; i++) {
//                                for (int j = 0; j < nBasis; j++) {
//                                    if (eqn.isDiffusive()) {
//                                        ADiag[nBasis * m + i][nBasis * q + j] -= Jac * weight * derfvR * baseLeft[i] * baseLeft[j];
//                                        if (eqn.isIPFace(TT[k])) {
//                                            ADiag[nBasis * m + i][nBasis * q + j] += Jac * weight * derIP * baseLeft[i] * baseLeft[j];
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        V[m] = 0;
//                    }
//                    double[] dV = new double[nEqs * dim];
//                    for (int d = 0; d < dim; d++) {
//                        for (int m = 0; m < nEqs; m++) {
//                            dV[nEqs * d + m] = h;
//                            double[] dWLh = Mat.plusVec(dWL, dV);
//                            double[] dWch = new double[nEqs * dim];
//                            for (int q = 0; q < nEqs; q++) {
//
//                                for (int r = 0; r < dim; r++) {
//                                    dWch[nEqs * r + q] = (dWLh[nEqs * r + q] + dWR[nEqs * r + q]) / 2 + par.beta0 * (WR[q] - WL[q]) / dL * n[k][p][r];
//                                }
//                            }
//                            double[] fvRh = eqn.numericalDiffusiveFlux(Wc, dWch, n[k][p], TT[k], elemData);
//                            dV[nEqs * d + m] = 0;
//                            for (int q = 0; q < nEqs; q++) {
//                                double derfvR = (fvRh[q] - fvR0[q]) / h;
//                                for (int i = 0; i < nBasis; i++) {
//                                    for (int j = 0; j < nBasis; j++) {
//                                        if (eqn.isDiffusive()) {
//                                            ADiag[nBasis * m + i][nBasis * q + j] -= Jac * weight * derfvR * dBaseLeft[i][d] * baseLeft[j];
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    // tato funkce vypocitava reziduum__________________________________________
//    private void residuum(double[] V, double[] K, double[][] KR) {
//
//        // vypocet toku hranici
//        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
//            residuumWall(k, V, K, KR[k]);
//        }
//
//        if (elemType.order > 1) { // volume integral only for DGFEM
//            double[] nor = new double[dim];
//            double[][] f = null;
//            double[][] fv = null;
//            double[] product = null;
//
//            for (int p = 0; p < Int.nIntVolume; p++) {
//                double[] base = Int.basisVolume[p];
//                double[][] dBase = Int.dXbasisVolume[p];
//                double Jac = Int.JacobianVolume[p];
//                double weight = Int.weightsVolume[p];
//
//                // interpolation of mesh velocity, and other data
//                double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[p]);
//
//                double[] WInt = new double[nEqs];
//                double[] dWInt = new double[dim * nEqs];
//                for (int m = 0; m < nEqs; m++) {
//                    for (int j = 0; j < nBasis; j++) {
//                        WInt[m] += (W[m * nBasis + j] + V[m * nBasis + j]) * base[j];
//                        for (int d = 0; d < dim; d++) {
//                            dWInt[nEqs * d + m] += (W[m * nBasis + j] + V[m * nBasis + j]) * dBase[j][d];
//                        }
//                    }
//                }
//                // convection
//                if (eqn.isConvective()) {
//                    f = new double[dim][];
//                    for (int d = 0; d < dim; d++) {
//                        nor[d] = 1;
//                        f[d] = eqn.convectiveFlux(WInt, nor, elemData);
//                        nor[d] = 0;
//                    }
//                }
//                // diffusion
//                if (eqn.isDiffusive()) {
//                    fv = new double[dim][];
//                    for (int d = 0; d < dim; d++) {
//                        nor[d] = 1;
//                        fv[d] = eqn.diffusiveFlux(WInt, dWInt, nor, elemData);
//                        nor[d] = 0;
//                    }
//                }
//                // production
//                if (eqn.isSourcePresent()) {
//                    product = eqn.sourceTerm(WInt, dWInt, elemData);
//                }
//
//                for (int m = 0; m < nEqs; m++) {
//                    for (int j = 0; j < nBasis; j++) {
//                        if (eqn.isConvective()) {
//                            double fsum = 0;
//                            double dWsum = 0;
//                            for (int d = 0; d < dim; d++) {
//                                fsum += (f[d][m] - u[d] * WInt[m]) * dBase[j][d];
//                                dWsum += dWInt[nEqs * d + m] * dBase[j][d];
//                            }
//                            K[nBasis * m + j] += Jac * weight * fsum - (eps + par.dampConst + dampInner[m]) * Jac * weight * dWsum;
//                        }
//                        if (eqn.isDiffusive()) {
//                            double fvsum = 0;
//                            for (int d = 0; d < dim; d++) {
//                                fvsum += fv[d][m] * dBase[j][d];
//                            }
//                            K[nBasis * m + j] -= Jac * weight * fvsum;
//                        }
//                        if (eqn.isSourcePresent()) {
//                            K[nBasis * m + j] += Jac * weight * product[m] * base[j];
//                        }
//                    }
//                }
//            }
//        } else// production term for FVM
//        {
//            if (eqn.isSourcePresent()) {
//                double[] Jac = Int.JacobianVolume;
//                double[] weights = Int.weightsVolume;
//
//                // interpolation of mesh velocity
//                double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[0]);
//
//                double[] WInt = new double[nEqs];
//                for (int j = 0; j < nEqs; j++) {
//                    WInt[j] = W[j] + V[j];
//
//                }
//                double[] dWInt = volumeDerivative(V, u, elemData);
//
//                // production
//                double[] product = eqn.sourceTerm(WInt, dWInt, elemData);
//
//                for (int m = 0; m < nEqs; m++) {
//                    if (eqn.isSourcePresent()) {
//                        K[m] += Jac[0] * weights[0] * product[m];
//                    }
//                }
//            }
//        }
//    }
//
//    // tato funkce vypocitava reziduum__________________________________________
//    private void residuumWall(int k, double[] V, double[] K, double[] KR) {
//        if (KR != null) {
//            Arrays.fill(KR, 0.0);
//        }
//        double[] fn = null;
//        double[] fvn = null;
//        int[] edgeIndex = Int.faces[k].faceIndexes;
//
//        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
//            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
//            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
//            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
//            double Jac = Int.faces[k].JacobianFace[p];
//            double weight = Int.faces[k].weightsFace[p];
//            double[] baseRight = null;
//            double[][] dBaseRight = null;
//            if (TT[k] > -1) {
//                baseRight = Int.faces[k].basisFaceRight[p];
//                dBaseRight = Int.faces[k].dXbasisFaceRight[p];
//            }
//
//            double dL = 0;
//            for (int d = 0; d < dim; d++) {
//                if (TT[k] > -1) {
//                    if (elems[TT[k]].insideComputeDomain) {
//                        dL += (Xs[d] - elems[TT[k]].Xs[d]) * n[k][p][d];
//                    } else {
//                        dL += 2 * (Xs[d] - Xes[k][d]) * n[k][p][d];
//                    }
//                } else {
//                    dL += (Xs[d] - Xes[k][d]) * n[k][p][d];
//                }
//            }
//            dL = Math.abs(dL);
//
//            // interpolation of mesh velocity
//            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);
//
//            double[] WL = new double[nEqs];
//            double[] WR = new double[nEqs];
//            double[] dWL = new double[dim * nEqs];
//            double[] dWR = new double[dim * nEqs];
//
//            // values from boundary inlet (WL, dWL)
//            if (elemType.order > 1) { // Discontinuous Galerkin Method
//                for (int m = 0; m < nEqs; m++) {
//                    for (int j = 0; j < nBasis; j++) {
//                        WL[m] += (W[m * nBasis + j] + V[m * nBasis + j]) * baseLeft[j];
//                        for (int d = 0; d < dim; d++) {
//                            dWL[nEqs * d + m] += (W[m * nBasis + j] + V[m * nBasis + j]) * dBaseLeft[j][d];
//                        }
//                    }
//                }
//            } else { // Finite volume method
//                dWL = volumeDerivative(V, u, elemData);
//                double sigmaL = FVMlimiter(dWL, par.FVMlimiter);
//                for (int m = 0; m < nEqs; m++) {
//                    double dW = 0;
//                    for (int d = 0; d < dim; d++) {
//                        dW = dW + (Int.faces[k].coordinatesFace[p][d] - Xs[d]) * dWL[nEqs * d + m];
//                    }
//                    WL[m] = W[m] + V[m] + sigmaL * dW;
//                }
//            }
//
//            // values from boundary outlet (WR, dWR)
//            if (TT[k] > -1) {
//                if (elems[TT[k]].elemType.order > 1) { // Discontinuous Galerkin Method
//                    double[] WRp = elems[TT[k]].getW(par.isExplicit, tLTS);
//                    for (int m = 0; m < nEqs; m++) {
//                        int nRBasis = elems[TT[k]].nBasis;
//                        for (int j = 0; j < nRBasis; j++) {
//                            WR[m] += WRp[m * nRBasis + j] * baseRight[j];
//                            for (int d = 0; d < dim; d++) {
//                                dWR[nEqs * d + m] += WRp[m * nRBasis + j] * dBaseRight[j][d];
//                            }
//                        }
//                    }
//                } else { // Finite Volume Method
//                    if (elems[TT[k]].insideComputeDomain) {
//                        dWR = elems[TT[k]].volumeDerivative(null, u, elemData);
//                    } else {
//                        System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
//                    }
//                    if (elems[TT[k]].insideComputeDomain) {
//                        double[] WRp = elems[TT[k]].getW(par.isExplicit, tLTS);
//                        double sigmaR = elems[TT[k]].FVMlimiter(dWR, par.FVMlimiter);
//                        for (int j = 0; j < nEqs; j++) {
//                            double dW = 0;
//                            for (int d = 0; d < dim; d++) {
//                                dW = dW + (Int.faces[k].coordinatesFace[p][d] - elems[TT[k]].Xs[d]) * dWR[nEqs * d + j];
//                            }
//                            WR[j] = WRp[j] + sigmaR * dW;
//                        }
//                    } else {
//                        System.arraycopy(elems[TT[k]].W, 0, WR, 0, nEqs);
//                    }
//                }
//            } else {
//                WR = eqn.boundaryValue(WL, n[k][p], TT[k], elemData);
//                System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
//            }
//
//            // inviscid flux in integration point
//            double vn = 0;
//            double[] Wale = new double[nEqs];
//            if (eqn.isConvective()) {
//                for (int d = 0; d < dim; d++) {
//                    vn = vn + u[d] * n[k][p][d];
//                }
//                if (TT[k] > -1) {
//                    for (int j = 0; j < nEqs; j++) {
//                        Wale[j] = (WL[j] + WR[j]) / 2;
//                    }
//                } else {
//                    System.arraycopy(WR, 0, Wale, 0, nEqs);
//                }
//                fn = eqn.numericalConvectiveFlux(WL, WR, n[k][p], TT[k], elemData);	// nevazky tok
//            }
//
//            // viscid flux in integration point
//            // DDG
//            if (eqn.isDiffusive()) {
//                // DDG
//                double beta0 = par.beta0;
//                if (par.order == 1) {
//                    beta0 = 0;
//                }
//                double[] Wc = new double[nEqs];
//                double[] dWc = new double[nEqs * dim];
//                for (int m = 0; m < nEqs; m++) {
//                    if (TT[k] > -1) {
//                        Wc[m] = (WL[m] + WR[m]) / 2;
//                    } else {
//                        Wc[m] = WR[m];
//                    }
//                    for (int d = 0; d < dim; d++) {
//                        dWc[nEqs * d + m] = (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2 + beta0 * (WR[m] - WL[m]) / dL * n[k][p][d];
//                    }
//                }
//                fvn = eqn.numericalDiffusiveFlux(Wc, dWc, n[k][p], TT[k], elemData);
//            }
//
//            for (int m = 0; m < nEqs; m++) {
//                double dWsum = 0;
//                if (TT[k] > -1) {
//                    for (int d = 0; d < dim; d++) {
//                        dWsum += (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2 * n[k][p][d];
//                    }
//                }
//                for (int j = 0; j < nBasis; j++) {
//                    double jwb = Jac * weight * baseLeft[j];
//                    if (eqn.isConvective()) {
//                        K[nBasis * m + j] -= jwb * (fn[m] - vn * Wale[m]);
//                        if (TT[k] > -1) {
//                            K[nBasis * m + j] += (0.5 * (eps + elems[TT[k]].eps) + par.dampConst) * jwb * dWsum;
//                        }
//                    }
//                    if (eqn.isDiffusive()) {
//                        K[nBasis * m + j] += jwb * fvn[m];
//                        if (TT[k] > -1) {
//                            K[nBasis * m + j] -= c_IP * jwb * (WL[m] - WR[m]);
//                        } else if (eqn.isIPFace(TT[k])) {
//                            K[nBasis * m + j] -= c_IP * jwb * (WL[m] - WR[m]);
//                        }
//                    }
//                }
//                if (KR != null) {
//                    int nRBasis = elems[TT[k]].nBasis;
//                    for (int j = 0; j < nRBasis; j++) {
//                        double jwb = Jac * weight * baseRight[j];
//                        if (eqn.isConvective()) {
//                            KR[nRBasis * m + j] -= jwb * (fn[m] - vn * Wale[m]);
//                            KR[nRBasis * m + j] += (0.5 * (eps + elems[TT[k]].eps) + par.dampConst) * jwb * dWsum;
//                        }
//                        if (eqn.isDiffusive()) {
//                            KR[nRBasis * m + j] += jwb * fvn[m];
//                            KR[nRBasis * m + j] -= c_IP * jwb * (WL[m] - WR[m]);
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    double[] volumeDerivative(double[] V, double[] u, ElementData elemData) {
//        double[] dW = new double[dim * nEqs];
//        double[] WL = new double[nEqs];
//        double[] WWall;
//        if (V != null) {
//            for (int j = 0; j < nEqs; j++) {
//                WL[j] = W[j] + V[j];
//            }
//        } else {
//            System.arraycopy(W, 0, WL, 0, nEqs);
//
//        }
//        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
//            double[] WR = new double[nEqs];
//            if (TT[k] > -1) {
//                if (elems[TT[k]].elemType.order > 1) { // DGFEM neigbhour
//                    if (elems[TT[k]].insideComputeDomain) {
//                        WR = elems[TT[k]].calculateAvgW();
//                    } else {
//                        System.arraycopy(WL, 0, WR, 0, nEqs);
//                    }
//                } else { // FVM neigbhour
//                    WR = elems[TT[k]].getW(par.isExplicit, tLTS);
//                }
//                WWall = Mat.times(Mat.plusVec(WR, WL), 0.5);
//            } else {
//                WWall = eqn.boundaryValue(WL, n[k][0], TT[k], elemData);
//            }
//
//            for (int j = 0; j < nEqs; j++) {
//                for (int d = 0; d < dim; d++) {
//                    dW[nEqs * d + j] += S[k] * WWall[j] * n[k][0][d] / area;
//                }
//            }
//        }
//
//        return dW;
//    }
//
//    double FVMlimiter(double[] dW, String limiterType) {
//        double sigma = 0;
//        double[] Wmin;
//        double[] Wmax;
//        switch (limiterType) {
//            case "noLimiter":
//                sigma = 1;
//                break;
//            case "barth":
//                Wmin = new double[nEqs];
//                Wmax = new double[nEqs];
//                for (int j = 0; j < nEqs; j++) {
//                    Wmin[j] = 1e7;
//                    Wmax[j] = -1e7;
//                }
//                for (int k = 0; k < nFaces; k++) {
//                    if (TT[k] > -1) {
//                        if (elems[TT[k]].insideComputeDomain) {
//                            double[] WR = elems[TT[k]].calculateAvgW();
//                            for (int j = 0; j < nEqs; j++) {
//                                if (WR[j] > Wmax[j]) {
//                                    Wmax[j] = WR[j];
//                                }
//                                if (WR[j] < Wmin[j]) {
//                                    Wmin[j] = WR[j];
//                                }
//                            }
//                        } else {
//                            return 0;
//                        }
//                    }
//                }
//                sigma = 2;
//                for (int k = 0; k < nFaces; k++) {
//                    for (int j = 0; j < nEqs; j++) {
//                        double a = 1;
//                        double deltaW = 0;
//                        for (int d = 0; d < dim; d++) {
//                            deltaW = deltaW + (Xes[k][d] - Xs[d]) * dW[nEqs * d + j];
//                        }
//                        double WL = W[j] + deltaW;
//                        if (WL > W[j]) {
//                            a = (Wmax[j] - W[j]) / (WL - W[j]);
//                        } else if (WL < W[j]) {
//                            a = (Wmin[j] - W[j]) / (WL - W[j]);
//                        }
//                        if (a < 0) {
//                            return 0;
//                        }
//                        if (a < sigma) {
//                            sigma = a;
//                        }
//                    }
//                }
//                break;
//
//            case "venka": // Venkatakrishnan
//                Wmin = new double[nEqs];
//                Wmax = new double[nEqs];
//                for (int j = 0; j < nEqs; j++) {
//                    Wmin[j] = 1e7;
//                    Wmax[j] = -1e7;
//                }
//                for (int k = 0; k < nFaces; k++) {
//                    if (TT[k] > -1) {
//                        if (elems[TT[k]].insideComputeDomain) {
//                            double[] WR = elems[TT[k]].calculateAvgW();
//                            for (int j = 0; j < nEqs; j++) {
//                                if (WR[j] > Wmax[j]) {
//                                    Wmax[j] = WR[j];
//                                }
//                                if (WR[j] < Wmin[j]) {
//                                    Wmin[j] = WR[j];
//                                }
//                            }
//                        } else {
//                            return 0;
//                        }
//                    }
//                }
//                sigma = 1;
//                for (int k = 0; k < nFaces; k++) {
//                    for (int j = 0; j < nEqs; j++) {
//                        double a = 1;
//                        double deltaW = 0;
//                        for (int d = 0; d < dim; d++) {
//                            deltaW = deltaW + (Xes[k][d] - Xs[d]) * dW[nEqs * d + j];
//                        }
//                        double WL = W[j] + deltaW;
//                        if (WL > W[j]) {
//                            a = venkatakrishnan((Wmax[j] - W[j]) / (WL - W[j]));
//                        } else if (WL < W[j]) {
//                            a = venkatakrishnan((Wmin[j] - W[j]) / (WL - W[j]));
//                        }
//                        if (a < 0) {
//                            return 0;
//                        }
//                        if (a < sigma) {
//                            sigma = a;
//                        }
//                    }
//                }
//                break;
//        }
//        return sigma;
//    }
//
//    double venkatakrishnan(double y) {
//        return (y * y + 2 * y) / (y * y + y + 2);
//    }
//
//    double[] interpolateVelocityAndFillElementDataObjectOnVolume(double[] innerInterpolant) {
//        // interpolation of mesh velocity, and other data
//        double[] u = new double[dim];
//        double[] currentXi = new double[dim];
//        elemData.currentWallDistance = 0;
//        for (int j = 0; j < nVertices; j++) {
//            for (int d = 0; d < dim; d++) {
//                u[d] += innerInterpolant[j] * U[j][d];
//                //elemData.currentX[d] += innerInterpolant[j] * vertices[j][d];
//                currentXi[d] += innerInterpolant[j] * transform.coordsXi[j][d];
//            }
//            elemData.currentWallDistance += innerInterpolant[j] * wallDistance[j];
//        }
//        elemData.currentX = transform.getX(currentXi);
//
//        if (par.externalField) {
//            int nExternalField = externalField[0].length;
//            elemData.externalField = new double[nExternalField];
//            for (int j = 0; j < nVertices; j++) {
//                for (int r = 0; r < nExternalField; r++) {
//                    elemData.externalField[r] += innerInterpolant[j] * externalField[j][r];
//                }
//            }
//        }
//        elemData.currentT = mesh.t;
//        elemData.integralMonitor = mesh.integralMonitor;
//        elemData.elemIndex = index;
//        if (par.solutionAverage) {
//            elemData.Wavg = calculateAvgW();
//        }
//
//        return u;
//    }
//
//    double[] interpolateVelocityAndFillElementDataObjectOnFace(int k, double[] innerInterpolant, int[] edgeIndex) {
//        // interpolation of mesh velocity
//        double[] u = new double[dim];
//        double[] currentXi = new double[dim];
//        elemData.currentWallDistance = 0;
//        for (int j = 0; j < Int.faces[k].nVerticesEdge; j++) {
//            for (int d = 0; d < dim; d++) {
//                u[d] += innerInterpolant[j] * U[edgeIndex[j]][d];
//                elemData.meshVelocity[d] = u[d];
//                //elemData.currentX[d] += innerInterpolant[j] * vertices[edgeIndex[j]][d];
//                currentXi[d] += innerInterpolant[j] * transform.coordsXi[edgeIndex[j]][d];
//            }
//            elemData.currentWallDistance += innerInterpolant[j] * wallDistance[edgeIndex[j]];
//        }
//        elemData.currentX = transform.getX(currentXi);
//
//        if (par.externalField) {
//            int nExternalField = externalField[0].length;
//            elemData.externalField = new double[nExternalField];
//            for (int j = 0; j < Int.faces[k].nVerticesEdge; j++) {
//                for (int r = 0; r < nExternalField; r++) {
//                    elemData.externalField[r] += innerInterpolant[j] * externalField[edgeIndex[j]][r];
//                }
//            }
//        }
//        elemData.currentT = mesh.t;
//        elemData.integralMonitor = mesh.integralMonitor;
//        elemData.elemIndex = index;
//        if (par.solutionAverage) {
//            elemData.Wavg = calculateAvgW();
//        }
//
//        return u;
//    }
//
//    public void updateRHS(double[] x
//    ) {
//        int n = nEqs * nBasis;
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                RHS_loc[i] = RHS_loc[i] - ADiag[j][i] * x[gi_U[j]];
//            }
//            for (int k = 0; k < nFaces; k++) {
//                if (TT[k] > -1) {
//                    for (int j = 0; j < nEqs * elems[TT[k]].nBasis; j++) {
//                        RHS_loc[i] = RHS_loc[i] - ANeighs[k].A[j][i] * x[elems[TT[k]].gi_U[j]];
//                    }
//                }
//            }
//        }
//    }
//
//    /**
//     *
//     * @param dt
//     * @return L1norm(W - Wo)
//     */
//    public double calculateResiduumW(double dt
//    ) {
//        double rez = 0;
//        for (int m = 0; m < nEqs; m++) {
//            for (int j = 0; j < nBasis; j++) {
//                rez = rez + Math.abs(W[m * nBasis + j] - Wo[m * nBasis + j]) / dt;
//            }
//        }
//
//        return rez;
//    }
//
//    // ulozeni W do Wo
//    public void copyW2Wo() {
//        System.arraycopy(Wo, 0, Wo2, 0, nBasis * nEqs);
//        System.arraycopy(W, 0, Wo, 0, nBasis * nEqs);
//    }
//
//    // ulozeni Wo do W (pri potizich s resenim)
//    public void copyWo2W() {
//        System.arraycopy(Wo, 0, W, 0, nBasis * nEqs);
//    }
//
//    //__________________________________________________________________________
//    public double delta_t(double CFL) { //vypocet maximalniho casoveho kroku
//        double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(centreVolumeInterpolant);
//        double lam = eqn.maxEigenvalue(calculateAvgW(), elemData);
//        double dt = CFL * elemSize / lam;
//
//        return dt;
//    }
//
//    //__________________________________________________________________________
//    public double shock_senzor(double kap) {
//        double Se = 0;
//        double pod = 0;
//
//        double[][] base = Int.basisVolume;
//        double[] Jac = Int.JacobianVolume;
//        double[] weights = Int.weightsVolume;
//
//        double[] rhoTrunCoef = new double[nBasis];
//        if (elemType.order > 2) {
//            for (int i = 0; i < nBasis; i++) {
//                for (int j = 0; j < nBasis; j++) {
//                    rhoTrunCoef[i] += TrunOrd[i][j] * W[j];
//                }
//            }
//        } else {
//            double[] Ws = calculateAvgW();
//            if ("taylor".equals(basis.basisType)) {
//                rhoTrunCoef[0] = Ws[0];
//            } else {
//                Arrays.fill(rhoTrunCoef, Ws[0]);
//            }
//        }
//        for (int p = 0; p < Int.nIntVolume; p++) { // hodnoty funkci f a g v integracnich bodech
//            double rhoInt = 0;
//            double rhoTrun = 0;
//            for (int m = 0; m < nBasis; m++) {
//                rhoInt += W[m] * base[p][m];
//                rhoTrun += rhoTrunCoef[m] * base[p][m];
//            }
//            Se += Jac[p] * weights[p] * (rhoInt - rhoTrun) * (rhoInt - rhoTrun);
//            pod += Jac[p] * weights[p] * rhoInt * rhoInt;
//        }
//        Se = Math.log10(Math.abs(Se / pod));
//        double S0 = Math.log10(1. / Math.pow(elemType.order - 1, 4.0));
//        double shock = 0.5 * (1 + Math.sin(Math.PI * (Se - S0) / (2 * kap)));
//
//        if (Se < S0 - kap) {
//            shock = 0;
//        }
//        if (Se > S0 + kap) {
//            shock = 1;
//        }
//
//        if (Double.isNaN(shock)) {
//            return 0;
//        } else {
//            return shock;
//        }
//    }
//
//    //__________________________________________________________________________
//    public double[] shock_senzor_all_eqn(double kap) {
//
//        double[][] base = Int.basisVolume;
//        double[] Jac = Int.JacobianVolume;
//        double[] weights = Int.weightsVolume;
//
//        double[] shock = new double[nEqs];
//        for (int m = 0; m < nEqs; m++) {
//            double Se = 0;
//            double pod = 0;
//            double[] varTrunCoef = new double[nBasis];
//            if (elemType.order > 2) {
//                for (int i = 0; i < nBasis; i++) {
//                    for (int j = 0; j < nBasis; j++) {
//                        varTrunCoef[i] += TrunOrd[i][j] * W[m * nBasis + j];
//                    }
//                }
//            } else {
//                double[] Ws = calculateAvgW();
//                if ("taylor".equals(basis.basisType)) {
//                    varTrunCoef[0] = Ws[m];
//                } else {
//                    Arrays.fill(varTrunCoef, Ws[m]);
//                }
//            }
//            for (int p = 0; p < Int.nIntVolume; p++) { // hodnoty funkci f a g v integracnich bodech
//                double varInt = 0;
//                double varTrun = 0;
//                for (int i = 0; i < nBasis; i++) {
//                    varInt += W[m * nBasis + i] * base[p][i];
//                    varTrun += varTrunCoef[i] * base[p][i];
//                }
//                Se += Jac[p] * weights[p] * (varInt - varTrun) * (varInt - varTrun);
//                pod += Jac[p] * weights[p] * varInt * varInt;
//            }
//            Se = Math.log10(Math.abs(Se / pod));
//            double S0 = Math.log10(1. / Math.pow(elemType.order - 1, 4.0));
//            shock[m] = 0.5 * (1 + Math.sin(Math.PI * (Se - S0) / (2 * kap)));
//
//            if (Se < S0 - kap) {
//                shock[m] = 0;
//            }
//            if (Se > S0 + kap) {
//                shock[m] = 1;
//            }
//            if (Double.isNaN(shock[m])) {
//                shock[m] = 0;
//            }
//        }
//
//        return shock;
//    }
//
//    public double shock_senzor2() {
//        double shock = 0;
//        double sumS = 0;
//        double rhoLmax = 0;
//        for (int k = 0; k < nFaces; k++) {
//            for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
//                double[] baseLeft = Int.faces[k].basisFaceLeft[p];
//                double Jac = Int.faces[k].JacobianFace[p];
//                double weight = Int.faces[k].weightsFace[p];
//                double[] baseRight = null;
//                if (TT[k] > -1) {
//                    baseRight = Int.faces[k].basisFaceRight[p];
//                }
//
//                double rhoL = 0;
//                for (int j = 0; j < nBasis; j++) {
//                    rhoL += W[j] * baseLeft[j];
//                }
//                double rhoR = rhoL;
//                if (TT[k] > -1) {
//                    double[] WRp = elems[TT[k]].getW(par.isExplicit, tLTS);
//                    rhoR = 0;
//                    for (int j = 0; j < elems[TT[k]].nBasis; j++) {
//                        rhoL += WRp[j] * baseRight[j];
//                    }
//                }
//
//                shock += Jac * weight * (rhoL - rhoR);
//
//                if (rhoL > rhoLmax) {
//                    rhoLmax = rhoL;
//                }
//            }
//            sumS += S[k];
//        }
//        return shock / rhoLmax / sumS / Math.pow(elemSize, par.order / 2.0);
//    }
//
//    public int nastav_globalni_index_U(int s) {
//        gi_U = new int[nEqs * nBasis];
//        for (int i = 0; i < nBasis * nEqs; i++) {
//            gi_U[i] = s;
//            s = s + 1;
//        }
//        return s;
//    }
//
//    public void alocateNeigbourhsLinearSolver() {
//        ANeighs = new Neighbour[nFaces];
//        for (int k = 0; k < nFaces; k++) {
//            if (TT[k] > -1) {
//                ANeighs[k] = new Neighbour(TT[k], nBasis, elems[TT[k]].nBasis, nEqs, !par.isExplicit);
//            } else {
//                ANeighs[k] = new Neighbour(TT[k], nBasis, 0, nEqs, !par.isExplicit);
//            }
//        }
//    }
//
//    public double sqr() {
//        double nrm = 0;
//        for (int i = 0; i < nEqs * nBasis; i++) {
//            nrm = nrm + RHS_loc[i] * RHS_loc[i];
//        }
//        return nrm;
//    }
//
//    /**
//     * W = W + x
//     *
//     * @param x
//     */
//    public void updateW(double[] x) {
//        for (int i = 0; i < nEqs * nBasis; i++) {
//            W[i] += x[gi_U[i]];
//        }
//    }
//
//    // vypocet geometrie _______________________________________________________
//    public void computeGeometry() {
//        Xes = new double[nFaces][];
//        for (int k = 0; k < nFaces; k++) {
//            Xes[k] = Int.faces[k].faceTransform.getXs();
//        }
//        S = new double[nFaces];
//        for (int k = 0; k < nFaces; k++) {
//            for (int i = 0; i < Int.faces[k].nIntEdge; i++) {
//                S[k] = S[k] + Int.faces[k].JacobianFace[i] * Int.faces[k].weightsFace[i];
//            }
//        }
//
//        n = new double[nFaces][][];
//        for (int k = 0; k < nFaces; k++) {
//            n[k] = new double[Int.faces[k].nIntEdge][dim];
//            double[][] faceCoords = Int.faces[k].quadFace.getCoords();
//            for (int p = 0; p < Int.faces[k].nIntEdge; p++) {
//                n[k][p] = Int.faces[k].faceTransform.getNormal(faceCoords[p]);
//            }
//            //Mat.print(n[k][0]);
//        }
//
//        if (dim == 1) {  // only for 1D case
//            n[0][0][0] = -1;
//        }
//
//        Xs = transform.getXs();
//
//        area = 0;
//        for (int i = 0; i < Int.nIntVolume; i++) {
//            area = area + Int.JacobianVolume[i] * Int.weightsVolume[i];
//        }
//        //System.out.println(area);
//
//        elemSize = 0;
//        for (int k = 0; k < nFaces; k++) {
//            elemSize += S[k];
//        }
//        elemSize = nFaces * area / elemSize;
//
//        faceIndexReverse = new int[nFaces];
//        for (int k = 0; k < nFaces; k++) {
//            if (TT[k] > -1) {
//                for (int s = 0; s < elems[TT[k]].nFaces; s++) {
//                    if (elems[TT[k]].TT[s] == index) {
//                        faceIndexReverse[k] = s;
//                        break;
//                    }
//                }
//            }
//        }
//    }
//
//    /**
//     * Generovani matice hmotnosti.
//     *
//     * @param basis
//     */
//    public void constructMassMatrix() { // funkce pro generovani matic
//        // integracni vzorec pro vypocet matice hmotnosti musi mit prislusny rad, zkontrolovat!!!!!!!!!!   
//
//        double[][] base = Int.basisVolume;
//        double[] Jac = Int.JacobianVolume;
//        double[] weights = Int.weightsVolume;
//
//        // matice hmotnosti
//        M = new double[nBasis][nBasis];
//        Is = new double[nBasis];
//        for (int i = 0; i < nBasis; i++) {
//            for (int j = 0; j < nBasis; j++) {
//                for (int p = 0; p < Int.nIntVolume; p++) {
//                    M[i][j] = M[i][j] + Jac[p] * weights[p] * base[p][i] * base[p][j];
//                }
//            }
//            for (int p = 0; p < Int.nIntVolume; p++) {
//                Is[i] = Is[i] + Jac[p] * weights[p] * base[p][i] / area;
//            }
//        }
//
//        // initialization of mass matrixes
//        Mo = new double[nBasis][nBasis];
//        Mo2 = new double[nBasis][nBasis];
//        for (int i = 0; i < nBasis; i++) {
//            for (int j = 0; j < nBasis; j++) {
//                Mo[i][j] = M[i][j];
//                Mo2[i][j] = M[i][j];
//            }
//        }
//    }
//
//    public void recalculateMassMatrix() {
//        double[][] base = Int.basisVolume;
//        double[] Jac = Int.JacobianVolume;
//        double[] weights = Int.weightsVolume;
//
//        for (int i = 0; i < nBasis; i++) {
//            for (int j = 0; j < nBasis; j++) {
//                M[i][j] = 0;
//                for (int p = 0; p < Int.nIntVolume; p++) {
//                    M[i][j] = M[i][j] + Jac[p] * weights[p] * base[p][i] * base[p][j];
//                }
//            }
//            Is[i] = 0;
//            for (int p = 0; p < Int.nIntVolume; p++) {
//                Is[i] = Is[i] + Jac[p] * weights[p] * base[p][i] / area;
//            }
//        }
//    }
//
//    public void nextTimeLevelMassMatrixes() {
//        for (int i = 0; i < nBasis; i++) {
//            for (int j = 0; j < nBasis; j++) {
//                Mo2[i][j] = Mo[i][j];
//                Mo[i][j] = M[i][j];
//            }
//        }
//    }
//
//    public void computeOrderTruncationMatrix() { // truncation to linear order matrix
//        if (elemType.order > 2) {
//            double[][] base = Int.basisVolume;
//            double[] Jac = Int.JacobianVolume;
//            double[] weights = Int.weightsVolume;
//            double[][] Xcoord = Int.transform.getX(Int.quadVolume.getCoords());
//
//            int nLinearBasis = 1 + dim;
//
//            double[][] Mhat = new double[nLinearBasis][nLinearBasis];
//            double[][] T = new double[nLinearBasis][nBasis];
//
//            for (int i = 0; i < nLinearBasis; i++) {
//                for (int j = 0; j < nLinearBasis; j++) {
//                    for (int p = 0; p < Int.nIntVolume; p++) {
//                        double bi = taylorBase(i, Xcoord[p]);
//                        double bj = taylorBase(j, Xcoord[p]);
//                        Mhat[i][j] += Jac[p] * weights[p] * bi * bj;
//                    }
//                }
//                for (int j = 0; j < nBasis; j++) {
//                    for (int p = 0; p < Int.nIntVolume; p++) {
//                        double bi = taylorBase(i, Xcoord[p]);
//                        double bj = base[p][j];
//                        T[i][j] += Jac[p] * weights[p] * bi * bj;
//                    }
//                }
//            }
//            TrunOrd = Mat.times(Mat.invert(Mhat), T);
//            TrunOrd = Mat.times(Mat.transpose(T), TrunOrd);
//            TrunOrd = Mat.times(Mat.invert(M), TrunOrd);
//        }
//    }
//
//    double taylorBase(int i, double[] X) {
//        if (i > 0) {
//            return (X[i - 1] - Xs[i - 1]);
//        } else {
//            return 1;
//        }
//    }
//
//    //____________________________________________________________________________
//    public double[] calculateAvgW() { // funkce vraci stredni hodnotu W na kontrolnim objemu K
//        double[] Ws = new double[nEqs];
//        for (int j = 0; j < nEqs; j++) {
//            for (int m = 0; m < nBasis; m++) {
//                Ws[j] = Ws[j] + W[j * nBasis + m] * Is[m];
//            }
//        }
//        return Ws;
//    }
//
//    public Equation getEqs() {
//        return eqn;
//    }
//
//    public int getNEqs() {
//        return nEqs;
//    }
//
//    public double[] getW(boolean inTime, double tIn) {
//        if (inTime) { // LTS explicit
//            double[] Waux = new double[nBasis * nEqs];
//            double tInterp;
//            if (tLTS == tLTSold) {
//                tInterp = 1;
//            } else {
//                tInterp = (tIn - tLTSold) / (tLTS - tLTSold);
//            }
//            switch (stepLTS) {
//                case 0:
//                    for (int j = 0; j < nEqs * nBasis; j++) {
//                        Waux[j] = tInterp * W[j] + (1 - tInterp) * WLTS[j];
//                    }
//                    break;
//                case 1:
//                    for (int j = 0; j < nEqs * nBasis; j++) {
//                        Waux[j] = tInterp * W1LTS[j] + (1 - tInterp) * W1LTSo[j];
//                    }
//                    break;
//                case 2:
//                    for (int j = 0; j < nEqs * nBasis; j++) {
//                        Waux[j] = tInterp * W2LTS[j] + (1 - tInterp) * W2LTSo[j];
//                    }
//                    break;
//            }
//            return Waux;
//        } else { // implicit method
//            return W;
//        }
//    }
//
//    public boolean geometryCheck(boolean writeOut) {
//        boolean isOK = true;
//        double[] test1 = new double[dim];
//        double[] test2 = new double[dim];
//        double sumTest1 = 0;
//        double sumTest2 = 0;
//
//        for (int k = 0; k < nFaces; k++) {
//            Integration.Face face = Int.faces[k];
//            for (int p = 0; p < face.nIntEdge; p++) { // edge integral
//                double[] X = face.faceTransform.getX(face.quadFace.coords[p]);
//                for (int d = 0; d < dim; d++) {
//                    test1[d] += face.weightsFace[p] * face.JacobianFace[p] * n[k][p][d];
//                    test2[d] += X[d] * face.weightsFace[p] * face.JacobianFace[p] * n[k][p][d];
//                }
//            }
//        }
//        for (int d = 0; d < dim; d++) {
//            sumTest1 += test1[d];
//            sumTest2 += test2[d] / dim;
//            if (test2[d] < 0) {
//                if (writeOut) {
//                    System.out.println("Mesh element error, bad points order in file elements.txt!");
//                }
//                isOK = false;
//            }
//        }
//        double testTol2 = 1e-11;
//        if (elemType.order == 1) {
//            testTol2 = 1e-4;
//        }
//        if (Math.abs(sumTest1) > 1e-11 || Math.abs(sumTest2 - area) > testTol2 || sumTest2 < 0) {
//            if (writeOut) {
//                System.out.println("Mesh element error! Control sum 1: " + Math.abs(sumTest1) + ", control sum 2: " + Math.abs(sumTest2 - area));
//            }
//            isOK = false;
//        }
//
////            // check normals
////            boolean err = false;
////            for (int k = 0; k < nFaces; k++) {
////                Face face = Int.faces[k];
////                for (int p = 0; p < face.nIntEdge; p++) { // edge integral
////                    double[] X = face.faceTransform.getX(face.quadFace.coords[p]);
////                    double norOrient = Mat.scalar(n[k][p], Mat.minusVec(X, Xs));
////                    if (norOrient < 0) {
////                        err = true;
////
////                    }
////                }
////                if (err) {
////                    System.out.println("Mesh orientation error!");
////                }
////            }
//        return isOK;
//    }
//
//    // optimalisation
//    // for optimization toolbox, generate residuum(W)
//    public void exportLocalFunctionalDerivative() {
//        optimFunDer = new double[nBasis * nEqs];
//        double[] V = new double[nBasis * nEqs];
//        double Iw = optimalisationFunctional.combineFunctionals(computeFunctional(W));
//        for (int i = 0; i < nBasis * nEqs; i++) {
//            V[i] = par.h;
//            double Iwh = optimalisationFunctional.combineFunctionals(computeFunctional(Mat.plusVec(W, V)));
//            optimFunDer[i] = (Iwh - Iw) / par.h;
//            V[i] = 0;
//        }
//    }
//
//    // for optimization toolbox, generate residuum(W)
//    public void exportLocalR() {
//        double[] V = new double[nBasis * nEqs];
//        double[] Rw = new double[nBasis * nEqs];
//        double[][] RwN = new double[nFaces][];
//        residuum(V, Rw, RwN);
//        System.arraycopy(Rw, 0, RHS_loc, 0, Rw.length);
//    }
//
//    // for optimization toolbox, generate only residuum derivation (dR/dW)
//    public void exportLocalJacobiMatrix() {
//        // vnitrni element - krivkovy i objemovy integral
//        double[] V = new double[nBasis * nEqs];
//        double[] Rw = new double[nBasis * nEqs];
//        double[][] RwNeigh = new double[nFaces][];
//        double[][] RwNeighH = new double[nFaces][];
//        for (int k = 0; k < nFaces; k++) {
//            if (TT[k] > -1) {
//                RwNeigh[k] = new double[elems[TT[k]].nBasis * nEqs];
//                RwNeighH[k] = new double[elems[TT[k]].nBasis * nEqs];
//            }
//        }
//        // compute residuum
//        residuum(V, Rw, RwNeigh);
//
//        double h = par.h;
//        for (int i = 0; i < nBasis * nEqs; i++) {
//            for (int j = 0; j < Rw.length; j++) {
//                ADiag[i][j] = -Rw[j];
//            }
//            V[i] = h;
//            residuum(V, ADiag[i], RwNeighH);
//            V[i] = 0;
//            for (int k = 0; k < nFaces; k++) {
//                if (TT[k] > -1) {
//                    double[][] Aaux = elems[TT[k]].ANeighs[faceIndexReverse[k]].A;
//                    for (int j = 0; j < RwNeighH[k].length; j++) {
//                        Aaux[i][j] = (RwNeighH[k][j] - RwNeigh[k][j]) / h;
//                    }
//                }
//            }
//        }
//        Mat.divide(ADiag, -h);
//    }
//
//    public double[] computeFunctional(double[] V) {
//        if (V == null) {
//            V = new double[nBasis * nEqs];
//        }
//        int nF = optimalisationFunctional.getN();
//        double[] f = new double[nF];
//
//        // vypocet toku hranici
//        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
//            Mat.plusVecToVec(f, computeFunctionalsWall(k, V));
//        }
//
//        for (int p = 0; p < Int.nIntVolume; p++) {
//            double[] base = Int.basisVolume[p];
//            double[][] dBase = Int.dXbasisVolume[p];
//            double Jac = Int.JacobianVolume[p];
//            double weight = Int.weightsVolume[p];
//
//            // interpolation of mesh velocity
//            double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[p]);
//
//            double[] WInt = new double[nEqs];
//            double[] dWInt = new double[dim * nEqs];
//            for (int m = 0; m < nEqs; m++) {
//                for (int j = 0; j < nBasis; j++) {
//                    WInt[m] += (W[m * nBasis + j] + V[m * nBasis + j]) * base[j];
//                    for (int d = 0; d < dim; d++) {
//                        dWInt[nEqs * d + m] += (W[m * nBasis + j] + V[m * nBasis + j]) * dBase[j][d];
//                    }
//                }
//            }
//            double[] aux = optimalisationFunctional.insideValue(WInt, dWInt, elemData);
//            for (int i = 0; i < nF; i++) {
//                f[i] += Jac * weight * aux[i];
//            }
//        }
//        return f;
//    }
//
//    // tato funkce vypocitava reziduum__________________________________________
//    public double[] computeFunctionalsWall(int k, double[] V) {
//        int nF = optimalisationFunctional.getN();
//        double[] f = new double[nF];
//        int[] edgeIndex = Int.faces[k].faceIndexes;
//
//        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
//            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
//            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
//            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
//            double Jac = Int.faces[k].JacobianFace[p];
//            double weight = Int.faces[k].weightsFace[p];
//
//            // interpolation of mesh velocity
//            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);
//
//            double[] WL = new double[nEqs];
//            double[] dWL = new double[dim * nEqs];
//
//            // values from boundary inlet (WL, dWL)
//            for (int j = 0; j < nEqs; j++) {
//                for (int m = 0; m < nBasis; m++) {
//                    WL[j] += (W[j * nBasis + m] + V[j * nBasis + m]) * baseLeft[m];
//                    for (int d = 0; d < dim; d++) {
//                        dWL[nEqs * d + j] += (W[j * nBasis + m] + V[j * nBasis + m]) * dBaseLeft[m][d];
//                    }
//                }
//            }
//            double[] aux = optimalisationFunctional.boundaryValue(WL, dWL, n[k][p], TT[k], elemData);
//            for (int i = 0; i < nF; i++) {
//                f[i] += Jac * weight * aux[i];
//            }
//        }
//        return f;
//    }
//
//    public double[] computeIntegralSolutionMonitor(int nIntegralMonitor) {
//        double[] f = new double[nIntegralMonitor];
//
//        // vypocet toku hranici
//        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
//            computeIntegralSolutionMonitorOnWall(k, f);
//        }
//
//        for (int p = 0; p < Int.nIntVolume; p++) {
//            double[] base = Int.basisVolume[p];
//            double[][] dBase = Int.dXbasisVolume[p];
//            double Jac = Int.JacobianVolume[p];
//            double weight = Int.weightsVolume[p];
//
//            // interpolation of mesh velocity
//            double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[p]);
//
//            double[] WInt = new double[nEqs];
//            double[] dWInt = new double[dim * nEqs];
//            for (int m = 0; m < nEqs; m++) {
//                for (int j = 0; j < nBasis; j++) {
//                    WInt[m] += W[m * nBasis + j] * base[j];
//                    for (int d = 0; d < dim; d++) {
//                        dWInt[nEqs * d + m] += W[m * nBasis + j] * dBase[j][d];
//                    }
//                }
//            }
//            double[] aux = mesh.solMonitor.insideValue(WInt, dWInt, elemData);
//            if (aux != null) {
//                for (int r = 0; r < aux.length; r++) {
//                    f[r] += Jac * weight * aux[r];
//                }
//            }
//        }
//        return f;
//    }
//
//    // tato funkce vypocitava reziduum__________________________________________
//    public void computeIntegralSolutionMonitorOnWall(int k, double[] f) {
//        int[] edgeIndex = Int.faces[k].faceIndexes;
//
//        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
//            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
//            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
//            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
//            double Jac = Int.faces[k].JacobianFace[p];
//            double weight = Int.faces[k].weightsFace[p];
//
//            // interpolation of mesh velocity
//            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);
//
//            double[] WL = new double[nEqs];
//            double[] dWL = new double[dim * nEqs];
//
//            // values from boundary inlet (WL, dWL)
//            for (int j = 0; j < nEqs; j++) {
//                for (int m = 0; m < nBasis; m++) {
//                    WL[j] += W[j * nBasis + m] * baseLeft[m];
//                    for (int d = 0; d < dim; d++) {
//                        dWL[nEqs * d + j] += W[j * nBasis + m] * dBaseLeft[m][d];
//                    }
//                }
//            }
//            double[] aux = mesh.solMonitor.boundaryValue(WL, dWL, n[k][p], TT[k], elemData);
//            if (aux != null) {
//                for (int r = 0; r < aux.length; r++) {
//                    f[r] += Jac * weight * aux[r];
//                }
//            }
//        }
//    }
//}
