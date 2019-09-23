/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.ElementData;
import flowpro.api.Mat;
import flowpro.core.Mesh;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.elementType.ElementType;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class DGFEM extends Element {
    
    public DGFEM(int index, double[][] vertices, double[][] Uinit, double[] wallDistance, double[][] externalField, int[] TT, int[] TP, int[] TEale, int[] TEshift, double[][] shift, FaceCurvature fCurv, double[][] blendFun, double[] initW,
            Mesh mesh, ElementType elemType) throws IOException {
        super(index, vertices, Uinit, wallDistance, externalField, TT, TP, TEale, TEshift, shift, fCurv, blendFun, initW, mesh, elemType);
    }
    
    public void initMethod() {
        computeOrderTruncationMatrix();
    }
    
    public void initCondition() {
        vertices0 = new double[nVertices][dim];
        for (int d = 0; d < dim; d++) {
            for (int i = 0; i < nVertices; i++) {
                vertices0[i][d] = vertices[i][d];
                U[i][d] = Uinit[i][d];
            }
        }

        // fill the solution vector with initial condition
        W = new double[nBasis * nEqs];
        Wo = new double[nBasis * nEqs];
        Wo2 = new double[nBasis * nEqs];
        if (initW.length == nBasis * nEqs) {
            System.arraycopy(initW, 0, W, 0, nBasis * nEqs);
            System.arraycopy(initW, 0, Wo, 0, nBasis * nEqs);
            System.arraycopy(initW, 0, Wo2, 0, nBasis * nEqs);
        } else {
            switch (basis.basisType) {
                case "lagrange":
                    for (int i = 0; i < nBasis; i++) {
                        for (int j = 0; j < nEqs; j++) {
                            W[j * nBasis + i] = initW[j];
                            Wo[j * nBasis + i] = initW[j];
                            Wo2[j * nBasis + i] = initW[j];
                        }
                    }
                    break;
                case "orthogonal":
                case "taylor":
                    for (int j = 0; j < nEqs; j++) {
                        W[j * nBasis] = initW[j];
                        Wo[j * nBasis] = initW[j];
                        Wo2[j * nBasis] = initW[j];
                    }
                    break;
            }
        }
    }
    
    // tato funkce vypocitava reziduum__________________________________________
    public void residuum(double[] V, double[] K, double[][] KR) {

        // vypocet toku hranici
        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
            residuumWall(k, V, K, KR[k]);
        }

        if (elemType.order > 1) { // volume integral only for DGFEM
            double[] nor = new double[dim];
            double[][] f = null;
            double[][] fv = null;
            double[] product = null;

            for (int p = 0; p < Int.nIntVolume; p++) {
                double[] base = Int.basisVolume[p];
                double[][] dBase = Int.dXbasisVolume[p];
                double Jac = Int.JacobianVolume[p];
                double weight = Int.weightsVolume[p];

                // interpolation of mesh velocity, and other data
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
                // convection
                if (eqn.isConvective()) {
                    f = new double[dim][];
                    for (int d = 0; d < dim; d++) {
                        nor[d] = 1;
                        f[d] = eqn.convectiveFlux(WInt, nor, elemData);
                        nor[d] = 0;
                    }
                }
                // diffusion
                if (eqn.isDiffusive()) {
                    fv = new double[dim][];
                    for (int d = 0; d < dim; d++) {
                        nor[d] = 1;
                        fv[d] = eqn.diffusiveFlux(WInt, dWInt, nor, elemData);
                        nor[d] = 0;
                    }
                }
                // production
                if (eqn.isSourcePresent()) {
                    product = eqn.sourceTerm(WInt, dWInt, elemData);
                }

                for (int m = 0; m < nEqs; m++) {
                    for (int j = 0; j < nBasis; j++) {
                        if (eqn.isConvective()) {
                            double fsum = 0;
                            double dWsum = 0;
                            for (int d = 0; d < dim; d++) {
                                fsum += (f[d][m] - u[d] * WInt[m]) * dBase[j][d];
                                dWsum += dWInt[nEqs * d + m] * dBase[j][d];
                            }
                            K[nBasis * m + j] += Jac * weight * fsum - (eps + par.dampConst + dampInner[m]) * Jac * weight * dWsum;
                        }
                        if (eqn.isDiffusive()) {
                            double fvsum = 0;
                            for (int d = 0; d < dim; d++) {
                                fvsum += fv[d][m] * dBase[j][d];
                            }
                            K[nBasis * m + j] -= Jac * weight * fvsum;
                        }
                        if (eqn.isSourcePresent()) {
                            K[nBasis * m + j] += Jac * weight * product[m] * base[j];
                        }
                    }
                }
            }
        } else// production term for FVM
        {
            if (eqn.isSourcePresent()) {
                double[] Jac = Int.JacobianVolume;
                double[] weights = Int.weightsVolume;

                // interpolation of mesh velocity
                double[] u = interpolateVelocityAndFillElementDataObjectOnVolume(Int.interpolantVolume[0]);

                double[] WInt = new double[nEqs];
                for (int j = 0; j < nEqs; j++) {
                    WInt[j] = W[j] + V[j];

                }
                double[] dWInt = volumeDerivative(V, u, elemData);

                // production
                double[] product = eqn.sourceTerm(WInt, dWInt, elemData);

                for (int m = 0; m < nEqs; m++) {
                    if (eqn.isSourcePresent()) {
                        K[m] += Jac[0] * weights[0] * product[m];
                    }
                }
            }
        }
    }

    // tato funkce vypocitava reziduum__________________________________________
    public void residuumWall(int k, double[] V, double[] K, double[] KR) {
        if (KR != null) {
            Arrays.fill(KR, 0.0);
        }
        double[] fn = null;
        double[] fvn = null;
        int[] edgeIndex = Int.faces[k].faceIndexes;

        for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
            double[] innerInterpolant = Int.faces[k].interpolantFace[p];
            double[] baseLeft = Int.faces[k].basisFaceLeft[p];
            double[][] dBaseLeft = Int.faces[k].dXbasisFaceLeft[p];
            double Jac = Int.faces[k].JacobianFace[p];
            double weight = Int.faces[k].weightsFace[p];
            double[] baseRight = null;
            double[][] dBaseRight = null;
            if (TT[k] > -1) {
                baseRight = Int.faces[k].basisFaceRight[p];
                dBaseRight = Int.faces[k].dXbasisFaceRight[p];
            }

            double dL = 0;
            for (int d = 0; d < dim; d++) {
                if (TT[k] > -1) {
                    if (elems[TT[k]].insideComputeDomain) {
                        dL += (Xs[d] - elems[TT[k]].Xs[d]) * n[k][p][d];
                    } else {
                        dL += 2 * (Xs[d] - Xes[k][d]) * n[k][p][d];
                    }
                } else {
                    dL += (Xs[d] - Xes[k][d]) * n[k][p][d];
                }
            }
            dL = Math.abs(dL);

            // interpolation of mesh velocity
            double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);

            double[] WL = new double[nEqs];
            double[] WR = new double[nEqs];
            double[] dWL = new double[dim * nEqs];
            double[] dWR = new double[dim * nEqs];

            // values from boundary inlet (WL, dWL)
            if (elemType.order > 1) { // Discontinuous Galerkin Method
                for (int m = 0; m < nEqs; m++) {
                    for (int j = 0; j < nBasis; j++) {
                        WL[m] += (W[m * nBasis + j] + V[m * nBasis + j]) * baseLeft[j];
                        for (int d = 0; d < dim; d++) {
                            dWL[nEqs * d + m] += (W[m * nBasis + j] + V[m * nBasis + j]) * dBaseLeft[j][d];
                        }
                    }
                }
            } else { // Finite volume method
                dWL = volumeDerivative(V, u, elemData);
                double sigmaL = FVMlimiter(dWL, par.FVMlimiter);
                for (int m = 0; m < nEqs; m++) {
                    double dW = 0;
                    for (int d = 0; d < dim; d++) {
                        dW = dW + (Int.faces[k].coordinatesFace[p][d] - Xs[d]) * dWL[nEqs * d + m];
                    }
                    WL[m] = W[m] + V[m] + sigmaL * dW;
                }
            }

            // values from boundary outlet (WR, dWR)
            if (TT[k] > -1) {
                if (elems[TT[k]].elemType.order > 1) { // Discontinuous Galerkin Method
                    double[] WRp = elems[TT[k]].W;
                    for (int m = 0; m < nEqs; m++) {
                        int nRBasis = elems[TT[k]].nBasis;
                        for (int j = 0; j < nRBasis; j++) {
                            WR[m] += WRp[m * nRBasis + j] * baseRight[j];
                            for (int d = 0; d < dim; d++) {
                                dWR[nEqs * d + m] += WRp[m * nRBasis + j] * dBaseRight[j][d];
                            }
                        }
                    }
                } else { // Finite Volume Method
                    if (elems[TT[k]].insideComputeDomain) {
                        dWR = ((DGFEM)elems[TT[k]]).volumeDerivative(null, u, elemData);
                    } else {
                        System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
                    }
                    if (elems[TT[k]].insideComputeDomain) {
                        double[] WRp = elems[TT[k]].W;
                        double sigmaR = ((DGFEM)elems[TT[k]]).FVMlimiter(dWR, par.FVMlimiter);
                        for (int j = 0; j < nEqs; j++) {
                            double dW = 0;
                            for (int d = 0; d < dim; d++) {
                                dW = dW + (Int.faces[k].coordinatesFace[p][d] - elems[TT[k]].Xs[d]) * dWR[nEqs * d + j];
                            }
                            WR[j] = WRp[j] + sigmaR * dW;
                        }
                    } else {
                        System.arraycopy(elems[TT[k]].W, 0, WR, 0, nEqs);
                    }
                }
            } else {
                WR = eqn.boundaryValue(WL, n[k][p], TT[k], elemData);
                System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
            }

            // inviscid flux in integration point
            double vn = 0;
            double[] Wale = new double[nEqs];
            if (eqn.isConvective()) {
                for (int d = 0; d < dim; d++) {
                    vn = vn + u[d] * n[k][p][d];
                }
                if (TT[k] > -1) {
                    for (int j = 0; j < nEqs; j++) {
                        Wale[j] = (WL[j] + WR[j]) / 2;
                    }
                } else {
                    System.arraycopy(WR, 0, Wale, 0, nEqs);
                }
                fn = eqn.numericalConvectiveFlux(WL, WR, n[k][p], TT[k], elemData);	// nevazky tok
            }

            // viscid flux in integration point
            // DDG
            if (eqn.isDiffusive()) {
                // DDG
                double beta0 = par.beta0;
                if (par.order == 1) {
                    beta0 = 0;
                }
                double[] Wc = new double[nEqs];
                double[] dWc = new double[nEqs * dim];
                for (int m = 0; m < nEqs; m++) {
                    if (TT[k] > -1) {
                        Wc[m] = (WL[m] + WR[m]) / 2;
                    } else {
                        Wc[m] = WR[m];
                    }
                    for (int d = 0; d < dim; d++) {
                        dWc[nEqs * d + m] = (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2 + beta0 * (WR[m] - WL[m]) / dL * n[k][p][d];
                    }
                }
                fvn = eqn.numericalDiffusiveFlux(Wc, dWc, n[k][p], TT[k], elemData);
            }

            for (int m = 0; m < nEqs; m++) {
                double dWsum = 0;
                if (TT[k] > -1) {
                    for (int d = 0; d < dim; d++) {
                        dWsum += (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2 * n[k][p][d];
                    }
                }
                for (int j = 0; j < nBasis; j++) {
                    double jwb = Jac * weight * baseLeft[j];
                    if (eqn.isConvective()) {
                        K[nBasis * m + j] -= jwb * (fn[m] - vn * Wale[m]);
                        if (TT[k] > -1) {
                            K[nBasis * m + j] += (0.5 * (eps + elems[TT[k]].eps) + par.dampConst) * jwb * dWsum;
                        }
                    }
                    if (eqn.isDiffusive()) {
                        K[nBasis * m + j] += jwb * fvn[m];
                        if (TT[k] > -1) {
                            K[nBasis * m + j] -= c_IP * jwb * (WL[m] - WR[m]);
                        } else if (eqn.isIPFace(TT[k])) {
                            K[nBasis * m + j] -= c_IP * jwb * (WL[m] - WR[m]);
                        }
                    }
                }
                if (KR != null) {
                    int nRBasis = elems[TT[k]].nBasis;
                    for (int j = 0; j < nRBasis; j++) {
                        double jwb = Jac * weight * baseRight[j];
                        if (eqn.isConvective()) {
                            KR[nRBasis * m + j] -= jwb * (fn[m] - vn * Wale[m]);
                            KR[nRBasis * m + j] += (0.5 * (eps + elems[TT[k]].eps) + par.dampConst) * jwb * dWsum;
                        }
                        if (eqn.isDiffusive()) {
                            KR[nRBasis * m + j] += jwb * fvn[m];
                            KR[nRBasis * m + j] -= c_IP * jwb * (WL[m] - WR[m]);
                        }
                    }
                }
            }
        }
    }
    
    
    double[] volumeDerivative(double[] V, double[] u, ElementData elemData) {
        double[] dW = new double[dim * nEqs];
        double[] WL = new double[nEqs];
        double[] WWall;
        if (V != null) {
            for (int j = 0; j < nEqs; j++) {
                WL[j] = W[j] + V[j];
            }
        } else {
            System.arraycopy(W, 0, WL, 0, nEqs);

        }
        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
            double[] WR = new double[nEqs];
            if (TT[k] > -1) {
                if (elems[TT[k]].elemType.order > 1) { // DGFEM neigbhour
                    if (elems[TT[k]].insideComputeDomain) {
                        WR = elems[TT[k]].calculateAvgW();
                    } else {
                        System.arraycopy(WL, 0, WR, 0, nEqs);
                    }
                } else { // FVM neigbhour
                    WR = elems[TT[k]].W;
                }
                WWall = Mat.times(Mat.plusVec(WR, WL), 0.5);
            } else {
                WWall = eqn.boundaryValue(WL, n[k][0], TT[k], elemData);
            }

            for (int j = 0; j < nEqs; j++) {
                for (int d = 0; d < dim; d++) {
                    dW[nEqs * d + j] += S[k] * WWall[j] * n[k][0][d] / area;
                }
            }
        }

        return dW;
    }

    double FVMlimiter(double[] dW, String limiterType) {
        double sigma = 0;
        double[] Wmin;
        double[] Wmax;
        switch (limiterType) {
            case "noLimiter":
                sigma = 1;
                break;
            case "barth":
                Wmin = new double[nEqs];
                Wmax = new double[nEqs];
                for (int j = 0; j < nEqs; j++) {
                    Wmin[j] = 1e7;
                    Wmax[j] = -1e7;
                }
                for (int k = 0; k < nFaces; k++) {
                    if (TT[k] > -1) {
                        if (elems[TT[k]].insideComputeDomain) {
                            double[] WR = elems[TT[k]].calculateAvgW();
                            for (int j = 0; j < nEqs; j++) {
                                if (WR[j] > Wmax[j]) {
                                    Wmax[j] = WR[j];
                                }
                                if (WR[j] < Wmin[j]) {
                                    Wmin[j] = WR[j];
                                }
                            }
                        } else {
                            return 0;
                        }
                    }
                }
                sigma = 2;
                for (int k = 0; k < nFaces; k++) {
                    for (int j = 0; j < nEqs; j++) {
                        double a = 1;
                        double deltaW = 0;
                        for (int d = 0; d < dim; d++) {
                            deltaW = deltaW + (Xes[k][d] - Xs[d]) * dW[nEqs * d + j];
                        }
                        double WL = W[j] + deltaW;
                        if (WL > W[j]) {
                            a = (Wmax[j] - W[j]) / (WL - W[j]);
                        } else if (WL < W[j]) {
                            a = (Wmin[j] - W[j]) / (WL - W[j]);
                        }
                        if (a < 0) {
                            return 0;
                        }
                        if (a < sigma) {
                            sigma = a;
                        }
                    }
                }
                break;

            case "venka": // Venkatakrishnan
                Wmin = new double[nEqs];
                Wmax = new double[nEqs];
                for (int j = 0; j < nEqs; j++) {
                    Wmin[j] = 1e7;
                    Wmax[j] = -1e7;
                }
                for (int k = 0; k < nFaces; k++) {
                    if (TT[k] > -1) {
                        if (elems[TT[k]].insideComputeDomain) {
                            double[] WR = elems[TT[k]].calculateAvgW();
                            for (int j = 0; j < nEqs; j++) {
                                if (WR[j] > Wmax[j]) {
                                    Wmax[j] = WR[j];
                                }
                                if (WR[j] < Wmin[j]) {
                                    Wmin[j] = WR[j];
                                }
                            }
                        } else {
                            return 0;
                        }
                    }
                }
                sigma = 1;
                for (int k = 0; k < nFaces; k++) {
                    for (int j = 0; j < nEqs; j++) {
                        double a = 1;
                        double deltaW = 0;
                        for (int d = 0; d < dim; d++) {
                            deltaW = deltaW + (Xes[k][d] - Xs[d]) * dW[nEqs * d + j];
                        }
                        double WL = W[j] + deltaW;
                        if (WL > W[j]) {
                            a = venkatakrishnan((Wmax[j] - W[j]) / (WL - W[j]));
                        } else if (WL < W[j]) {
                            a = venkatakrishnan((Wmin[j] - W[j]) / (WL - W[j]));
                        }
                        if (a < 0) {
                            return 0;
                        }
                        if (a < sigma) {
                            sigma = a;
                        }
                    }
                }
                break;
        }
        return sigma;
    }

    double venkatakrishnan(double y) {
        return (y * y + 2 * y) / (y * y + y + 2);
    }
    
    // Aplikace limiteru
    public void limiter(boolean isFirstIter) {
        eps = 0;
        c_IP = 0;
        if (elemType.order > 1) {
            // max eigenvalue
            double lam = eqn.maxEigenvalue(calculateAvgW(), elemData);

            // artificial viscosity
            if (par.dampTol > 0) {
                double g_shock = shock_senzor(par.dampTol);
                eps = lam * elemSize / elemType.order * g_shock;
            }

            if (par.dampInnerTol > 0) {
                // inner damping based on artificial viscosity
                double[] g_shockInner = eqn.combineShockSensors(shock_senzor_all_eqn(par.dampInnerTol));
                if (isFirstIter && !par.continueComputation) {
                    for (int i = 0; i < g_shockInner.length; i++) {
                        // g_shockInner[i] = 1;
                    }
                }

                for (int m = 0; m < nEqs; m++) {
                    double indicator = Math.max(g_shockInner[m], innerIndicatorOld[m]);
                    innerIndicatorOld[m] = 0.1 * indicator;
                    if (par.dampInnerCoef != null) {
                        if (par.dampInnerCoef.length == 1) {
                            dampInner[m] = par.dampInnerCoef[0] * indicator;
                        } else {
                            dampInner[m] = par.dampInnerCoef[m] * indicator;
                        }
                    } else {
                        dampInner[m] = lam * elemSize / elemType.order * indicator;
                    }
                }
            }

            if (eqn.isDiffusive()) {
                c_IP = par.penalty; // * elemSize;
            }
        }
    }
    
    //__________________________________________________________________________
    public double shock_senzor(double kap) {
        double Se = 0;
        double pod = 0;

        double[][] base = Int.basisVolume;
        double[] Jac = Int.JacobianVolume;
        double[] weights = Int.weightsVolume;

        double[] rhoTrunCoef = new double[nBasis];
        if (elemType.order > 2) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    rhoTrunCoef[i] += TrunOrd[i][j] * W[j];
                }
            }
        } else {
            double[] Ws = calculateAvgW();
            if ("taylor".equals(basis.basisType)) {
                rhoTrunCoef[0] = Ws[0];
            } else {
                Arrays.fill(rhoTrunCoef, Ws[0]);
            }
        }
        for (int p = 0; p < Int.nIntVolume; p++) { // hodnoty funkci f a g v integracnich bodech
            double rhoInt = 0;
            double rhoTrun = 0;
            for (int m = 0; m < nBasis; m++) {
                rhoInt += W[m] * base[p][m];
                rhoTrun += rhoTrunCoef[m] * base[p][m];
            }
            Se += Jac[p] * weights[p] * (rhoInt - rhoTrun) * (rhoInt - rhoTrun);
            pod += Jac[p] * weights[p] * rhoInt * rhoInt;
        }
        Se = Math.log10(Math.abs(Se / pod));
        double S0 = Math.log10(1. / Math.pow(elemType.order - 1, 4.0));
        double shock = 0.5 * (1 + Math.sin(Math.PI * (Se - S0) / (2 * kap)));

        if (Se < S0 - kap) {
            shock = 0;
        }
        if (Se > S0 + kap) {
            shock = 1;
        }

        if (Double.isNaN(shock)) {
            return 0;
        } else {
            return shock;
        }
    }

    //__________________________________________________________________________
    public double[] shock_senzor_all_eqn(double kap) {

        double[][] base = Int.basisVolume;
        double[] Jac = Int.JacobianVolume;
        double[] weights = Int.weightsVolume;

        double[] shock = new double[nEqs];
        for (int m = 0; m < nEqs; m++) {
            double Se = 0;
            double pod = 0;
            double[] varTrunCoef = new double[nBasis];
            if (elemType.order > 2) {
                for (int i = 0; i < nBasis; i++) {
                    for (int j = 0; j < nBasis; j++) {
                        varTrunCoef[i] += TrunOrd[i][j] * W[m * nBasis + j];
                    }
                }
            } else {
                double[] Ws = calculateAvgW();
                if ("taylor".equals(basis.basisType)) {
                    varTrunCoef[0] = Ws[m];
                } else {
                    Arrays.fill(varTrunCoef, Ws[m]);
                }
            }
            for (int p = 0; p < Int.nIntVolume; p++) { // hodnoty funkci f a g v integracnich bodech
                double varInt = 0;
                double varTrun = 0;
                for (int i = 0; i < nBasis; i++) {
                    varInt += W[m * nBasis + i] * base[p][i];
                    varTrun += varTrunCoef[i] * base[p][i];
                }
                Se += Jac[p] * weights[p] * (varInt - varTrun) * (varInt - varTrun);
                pod += Jac[p] * weights[p] * varInt * varInt;
            }
            Se = Math.log10(Math.abs(Se / pod));
            double S0 = Math.log10(1. / Math.pow(elemType.order - 1, 4.0));
            shock[m] = 0.5 * (1 + Math.sin(Math.PI * (Se - S0) / (2 * kap)));

            if (Se < S0 - kap) {
                shock[m] = 0;
            }
            if (Se > S0 + kap) {
                shock[m] = 1;
            }
            if (Double.isNaN(shock[m])) {
                shock[m] = 0;
            }
        }

        return shock;
    }

    public double shock_senzor2() {
        double shock = 0;
        double sumS = 0;
        double rhoLmax = 0;
        for (int k = 0; k < nFaces; k++) {
            for (int p = 0; p < Int.faces[k].nIntEdge; p++) { // edge integral
                double[] baseLeft = Int.faces[k].basisFaceLeft[p];
                double Jac = Int.faces[k].JacobianFace[p];
                double weight = Int.faces[k].weightsFace[p];
                double[] baseRight = null;
                if (TT[k] > -1) {
                    baseRight = Int.faces[k].basisFaceRight[p];
                }

                double rhoL = 0;
                for (int j = 0; j < nBasis; j++) {
                    rhoL += W[j] * baseLeft[j];
                }
                double rhoR = rhoL;
                if (TT[k] > -1) {
                    double[] WRp = elems[TT[k]].W;
                    rhoR = 0;
                    for (int j = 0; j < elems[TT[k]].nBasis; j++) {
                        rhoL += WRp[j] * baseRight[j];
                    }
                }

                shock += Jac * weight * (rhoL - rhoR);

                if (rhoL > rhoLmax) {
                    rhoLmax = rhoL;
                }
            }
            sumS += S[k];
        }
        return shock / rhoLmax / sumS / Math.pow(elemSize, par.order / 2.0);
    }
    
    public void computeMassMatrixAndBasisWeights() { // funkce pro generovani matic
        // integracni vzorec pro vypocet matice hmotnosti musi mit prislusny rad, zkontrolovat!!!!!!!!!!   

        double[][] base = Int.basisVolume;
        double[] Jac = Int.JacobianVolume;
        double[] weights = Int.weightsVolume;

        // matice hmotnosti
        M = new double[nBasis][nBasis];
        Is = new double[nBasis];
        for (int i = 0; i < nBasis; i++) {
            for (int j = 0; j < nBasis; j++) {
                for (int p = 0; p < Int.nIntVolume; p++) {
                    M[i][j] = M[i][j] + Jac[p] * weights[p] * base[p][i] * base[p][j];
                }
            }
            for (int p = 0; p < Int.nIntVolume; p++) {
                Is[i] = Is[i] + Jac[p] * weights[p] * base[p][i] / area;
            }
        }

        // initialization of mass matrixes
        for (int i = 0; i < nBasis; i++) {
            for (int j = 0; j < nBasis; j++) {
                Mo[i][j] = M[i][j];
                Mo2[i][j] = M[i][j];
            }
        }
    }

    public void recalculateMassMatrix() {
        double[][] base = Int.basisVolume;
        double[] Jac = Int.JacobianVolume;
        double[] weights = Int.weightsVolume;

        for (int i = 0; i < nBasis; i++) {
            for (int j = 0; j < nBasis; j++) {
                M[i][j] = 0;
                for (int p = 0; p < Int.nIntVolume; p++) {
                    M[i][j] = M[i][j] + Jac[p] * weights[p] * base[p][i] * base[p][j];
                }
            }
            Is[i] = 0;
            for (int p = 0; p < Int.nIntVolume; p++) {
                Is[i] = Is[i] + Jac[p] * weights[p] * base[p][i] / area;
            }
        }
    }
    
    public void computeOrderTruncationMatrix() { // truncation to linear order matrix
        if (elemType.order > 2) {
            double[][] base = Int.basisVolume;
            double[] Jac = Int.JacobianVolume;
            double[] weights = Int.weightsVolume;
            double[][] Xcoord = Int.transform.getX(Int.quadVolume.getCoords());

            int nLinearBasis = 1 + dim;

            double[][] Mhat = new double[nLinearBasis][nLinearBasis];
            double[][] T = new double[nLinearBasis][nBasis];

            for (int i = 0; i < nLinearBasis; i++) {
                for (int j = 0; j < nLinearBasis; j++) {
                    for (int p = 0; p < Int.nIntVolume; p++) {
                        double bi = taylorBase(i, Xcoord[p]);
                        double bj = taylorBase(j, Xcoord[p]);
                        Mhat[i][j] += Jac[p] * weights[p] * bi * bj;
                    }
                }
                for (int j = 0; j < nBasis; j++) {
                    for (int p = 0; p < Int.nIntVolume; p++) {
                        double bi = taylorBase(i, Xcoord[p]);
                        double bj = base[p][j];
                        T[i][j] += Jac[p] * weights[p] * bi * bj;
                    }
                }
            }
            TrunOrd = Mat.times(Mat.invert(Mhat), T);
            TrunOrd = Mat.times(Mat.transpose(T), TrunOrd);
            TrunOrd = Mat.times(Mat.invert(M), TrunOrd);
        }
    }

    double taylorBase(int i, double[] X) {
        if (i > 0) {
            return (X[i - 1] - Xs[i - 1]);
        } else {
            return 1;
        }
    }
}
