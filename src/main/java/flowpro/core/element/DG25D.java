///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package flowpro.core.element;
//
//import flowpro.api.FlowProProperties;
//import flowpro.api.Mat;
//import flowpro.core.Integration;
//import flowpro.core.Mesh;
//import flowpro.core.curvedBoundary.FaceCurvature;
//import flowpro.core.elementType.ElementType;
//import java.io.IOException;
//import java.util.Arrays;
//
///**
// *
// * @author obublik
// */
//public class DG25D extends Element {
//
//    public double penalty;    // interior penalty constant
//    public double beta0;  // direct discontinuous constant 
//
//    // artificial damping
//    public double dampTol;   // tolerance pro pridavnou viskozitu
//    public double dampInnerTol; // tolerance pro pridavnou viskozitu uvnitr
//    public double[] dampInnerCoef; // koeficienty pro pridavnou viskozitu uvnitr
//    public double dampConst; // konstantni pridavna viskozita 
//
//    // integration in z dimension
//    IntegrationXYZ IntXYZ;
//
//    public DG25D() {
//
//    }
//
//    public void set(int index, double[][] vertices, double[][] Uinit, double[] wallDistance, double[][] externalField, int[] TT, int[] TP, int[] TEale, int[] TEshift, double[][] shift, FaceCurvature fCurv, double[][] blendFun, double[] initW,
//            Mesh mesh, ElementType elemType) throws IOException {
//        super.set(index, vertices, Uinit, wallDistance, externalField, TT, TP, TEale, TEshift, shift, fCurv, blendFun, initW, mesh, elemType);
//    }
//
//    public void initMethod(FlowProProperties props) throws IOException {
//        computeOrderTruncationMatrix();
//
//        // damping
//        dampTol = props.getDouble("dampTol");
//        dampConst = props.getDouble("dampConst");
//
//        // iterative solver setting
//        if (props.containsKey("dampInnerTol")) {
//            dampInnerTol = props.getDouble("dampInnerTol");
//        } else {
//            dampInnerTol = 0;
//        }
//        if (props.containsKey("dampInnerCoef")) {
//            dampInnerCoef = props.getDoubleArray("dampInnerCoef");
//        } else {
//            dampInnerCoef = null;
//        }
//
//        if (props.containsKey("penalty")) {
//            penalty = props.getDouble("penalty");
//        } else {
//            penalty = 0;
//        }
//
//        if (props.containsKey("beta0")) {
//            beta0 = props.getDouble("beta0");
//        } else {
//            beta0 = 2;
//        }
//
//        IntXYZ = new IntegrationXYZ(Int, 1, 4);
//    }
//
//    public void initCondition() {
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
//    // tato funkce vypocitava reziduum__________________________________________
//    @Override
//    public void residuum(double[] V, double[] K, double[][] KR) {
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
//            for (int p = 0; p < IntXYZ.nIntVolume; p++) {
//                double[] base = IntXYZ.basisVolume[p];
//                double[][] dBase = IntXYZ.dXbasisVolume[p];
//                double Jac = IntXYZ.JacobianVolume[p];
//                double weight = IntXYZ.weightsVolume[p];
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
//                                fsum += f[d][m] * dBase[j][d];
//                                dWsum += dWInt[nEqs * d + m] * dBase[j][d];
//                            }
//                            K[nBasis * m + j] += Jac * weight * fsum - (eps + dampConst + dampInner[m]) * Jac * weight * dWsum;
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
//        }
//    }
//
//    // tato funkce vypocitava reziduum__________________________________________
//    @Override
//    public void residuumWall(int k, double[] V, double[] K, double[] KR) {
//        if (KR != null) {
//            Arrays.fill(KR, 0.0);
//        }
//        double[] fn = null;
//        double[] fvn = null;
//
//        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
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
//            double[] WL = new double[nEqs];
//            double[] WR = new double[nEqs];
//            double[] dWL = new double[dim * nEqs];
//            double[] dWR = new double[dim * nEqs];
//
//            // values from boundary inlet (WL, dWL)
//            for (int m = 0; m < nEqs; m++) {
//                for (int j = 0; j < nBasis; j++) {
//                    WL[m] += (W[m * nBasis + j] + V[m * nBasis + j]) * baseLeft[j];
//                    for (int d = 0; d < dim; d++) {
//                        dWL[nEqs * d + m] += (W[m * nBasis + j] + V[m * nBasis + j]) * dBaseLeft[j][d];
//                    }
//                }
//            }
//
//            // values from boundary outlet (WR, dWR)
//            if (TT[k] > -1) {
//                double[] WRp = elems[TT[k]].W;
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
//                fn = eqn.numericalConvectiveFlux(WL, WR, n[k][p], TT[k], elemData);	// nevazky tok
//            }
//
//            // viscid flux in integration point
//            // DDG
//            if (eqn.isDiffusive()) {
//                // DDG
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
//                        dWc[nEqs * d + m] = (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2;
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
//                        K[nBasis * m + j] -= jwb * fn[m];
//                        if (TT[k] > -1) {
//                            K[nBasis * m + j] += (0.5 * (eps + elems[TT[k]].eps) + dampConst) * jwb * dWsum;
////                        } else{
////                            System.out.println(u[1]);
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
//                            KR[nRBasis * m + j] -= jwb * fn[m];
//                            KR[nRBasis * m + j] += (0.5 * (eps + elems[TT[k]].eps) + dampConst) * jwb * dWsum;
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
//    // Aplikace limiteru
//    public void limiter(boolean isFirstIter) {
//        eps = 0;
//        c_IP = 0;
//        if (elemType.order > 1) {
//            // max eigenvalue
//            double lam = eqn.maxEigenvalue(calculateAvgW(), elemData);
//
//            // artificial viscosity
//            if (dampTol > 0) {
//                double g_shock = shock_senzor(dampTol);
//                eps = lam * elemSize / elemType.order * g_shock;
//            }
//
//            if (dampInnerTol > 0) {
//                // inner damping based on artificial viscosity
//                double[] g_shockInner = eqn.combineShockSensors(shock_senzor_all_eqn(dampInnerTol));
//                if (isFirstIter && !par.continueComputation) {
//                    for (int i = 0; i < g_shockInner.length; i++) {
//                        // g_shockInner[i] = 1;
//                    }
//                }
//
//                for (int m = 0; m < nEqs; m++) {
//                    double indicator = Math.max(g_shockInner[m], innerIndicatorOld[m]);
//                    innerIndicatorOld[m] = 0.1 * indicator;
//                    if (dampInnerCoef != null) {
//                        if (dampInnerCoef.length == 1) {
//                            dampInner[m] = dampInnerCoef[0] * indicator;
//                        } else {
//                            dampInner[m] = dampInnerCoef[m] * indicator;
//                        }
//                    } else {
//                        dampInner[m] = lam * elemSize / elemType.order * indicator;
//                    }
//                }
//            }
//
//            if (eqn.isDiffusive()) {
//                c_IP = penalty; // * elemSize;
//            }
//        }
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
//    public void computeMassMatrixAndBasisWeights() { // funkce pro generovani matic
//        // integracni vzorec pro vypocet matice hmotnosti musi mit prislusny rad, zkontrolovat!!!!!!!!!!   
//
//        double[][] base = Int.basisVolume;
//        double[] Jac = Int.JacobianVolume;
//        double[] weights = Int.weightsVolume;
//
//        // matice hmotnosti
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
//    @Override
//    public boolean isJacobiMatrixAssembly() {
//        return false;
//    }
//
//    @Override
//    public void residuumJacobi(double[][] ADiag, Implicit.Neighbour[] Sous) {
//    }
//
//    @Override
//    public void residuumWallJacobi(int k, double[][] ADiag, Implicit.Neighbour Sous) {
//    }
//
//    class IntegrationXYZ {
//
//        int nBasisZ;
//        int nIntZ;
//        int nIntVolume;
//        double[] weightsVolume;
//        double[][] basisVolume;
//        double[][][] dXibasisVolume;
//        double[] JacobianVolume;
//
//        IntegrationXYZ(Integration Int, int nBasisZ, int nIntZ) {
//            this.nIntZ = nIntZ;
//            this.nIntVolume = Int.nIntVolume * nIntZ;
//            this.nBasisZ = nBasisZ;
//
//            weightsVolume = new double[nIntVolume];
//            basisVolume = new double[nIntVolume][nBasisZ];
//            dXibasisVolume = new double[nIntVolume][nBasisZ][3];
//            for (int i = 0; i < nIntZ; i++) {
//                for (int j = 0; j < Int.nIntVolume; j++) {
//                    weightsVolume[i*Int.nIntVolume + j] = Int.weightsVolume[j] / nIntZ;
//                    double zeta = (i+0.5)/nIntZ;
//                    for(int k = 0; k < Int.basisVolume[0].length; k++){
//                        for(int m = 0; m < nBasisZ; m++){
//                            basisVolume[i*Int.nIntVolume + j][k*nBasisZ + m] = Int.basisVolume[j][k] * basis_z(m, zeta);
//                            
//                            dXibasisVolume[i*Int.nIntVolume + j][k*nBasisZ + m][0] = Int.dXibasisVolume[j][k][0] * basis_z(m, zeta);
//                            dXibasisVolume[i*Int.nIntVolume + j][k*nBasisZ + m][1] = Int.dXibasisVolume[j][k][1] * basis_z(m, zeta);
//                            dXibasisVolume[i*Int.nIntVolume + j][k*nBasisZ + m][2] = Int.basisVolume[j][k] * DX_basis_z(m, zeta);
//                        }
//                    }
//                }
//            }
//        }
//
//        double basis_z(int n, double zeta) {
//            double b = 0;
//            switch (n) {
//                case (0):
//                    b = 1;
//                    break;
//                case (1):
//                    b = Math.sin(2 * Math.PI * zeta);
//                    break;
//                case (2):
//                    b = Math.cos(2 * Math.PI * zeta);
//                    break;
//                case (3):
//                    b = Math.sin(4 * Math.PI * zeta);
//                    break;
//                case (4):
//                    b = Math.cos(4 * Math.PI * zeta);
//                    break;
//            }
//            return b;
//        }
//        
//        double DX_basis_z(int n, double zeta) {
//            double b = 0;
//            switch (n) {
//                case (0):
//                    b = 1;
//                    break;
//                case (1):
//                    b = 2 * Math.PI * Math.cos(2 * Math.PI * zeta);
//                    break;
//                case (2):
//                    b = -2 * Math.PI * Math.sin(2 * Math.PI * zeta);
//                    break;
//                case (3):
//                    b = 4 * Math.PI * Math.cos(4 * Math.PI * zeta);
//                    break;
//                case (4):
//                    b = -4 * Math.PI * Math.sin(4 * Math.PI * zeta);
//                    break;
//            }
//            return b;
//        }
//    }
//}
