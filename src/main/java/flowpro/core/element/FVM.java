/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
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
public class FVM extends Element {

    // Finite volume method limiter
    public String FVMlimiter;

    public FVM(){
        maxMethodSpatialOrder = 1;
    }
    
    public void set(int index, double[][] vertices, double[][] Uinit, double[] wallDistance, double[][] externalField, int[] TT, int[] TP, int[] TEale, int[] TEshift, double[][] shift, FaceCurvature fCurv, double[][] blendFun, double[] initW,
            Mesh mesh, ElementType elemType) throws IOException {
        super.set(index, vertices, Uinit, wallDistance, externalField, TT, TP, TEale, TEshift, shift, fCurv, blendFun, initW, mesh, elemType);
    }

    public void initMethod(FlowProProperties props) throws IOException {
        if (props.containsKey("FVMlimiter")) {
            FVMlimiter = props.getString("FVMlimiter");
        } else {
            FVMlimiter = "none";
        }
    }

    public void initCondition() {
        // fill the solution vector with initial condition
        W = new double[nEqs];
        Wo = new double[nEqs];
        Wo2 = new double[nEqs];
        if (initW.length == nEqs) {
            System.arraycopy(initW, 0, W, 0, nEqs);
            System.arraycopy(initW, 0, Wo, 0, nEqs);
            System.arraycopy(initW, 0, Wo2, 0, nEqs);
        } else {
            for (int j = 0; j < nEqs; j++) {
                W[j] = initW[j];
                Wo[j] = initW[j];
                Wo2[j] = initW[j];
            }
        }
    }

    // tato funkce vypocitava reziduum__________________________________________
    public void residuum(double[] V, double[] K, double[][] KR) {

        // vypocet toku hranici
        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
            residuumWall(k, V, K, KR[k]);
        }

        if (eqn.isSourcePresent()) {
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
                    K[m] += area * product[m];
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
        double[] innerInterpolant = Int.faces[k].interpolantFace[0];

        // interpolation of mesh velocity
        double[] u = interpolateVelocityAndFillElementDataObjectOnFace(k, innerInterpolant, edgeIndex);

        double[] WL = new double[nEqs];
        double[] WR = new double[nEqs];
        double[] dWR = new double[dim * nEqs];

        // values from boundary inlet (WL, dWL)
        double[] dWL = volumeDerivative(V, u, elemData);
        double sigmaL = FVMlimiter(dWL, FVMlimiter);
        for (int m = 0; m < nEqs; m++) {
            double dW = 0;
            for (int d = 0; d < dim; d++) {
                dW = dW + (Int.faces[k].coordinatesFace[0][d] - Xs[d]) * dWL[nEqs * d + m];
            }
            WL[m] = W[m] + V[m] + sigmaL * dW;
        }

        // values from boundary outlet (WR, dWR)
        if (TT[k] > -1) {
            if (elems[TT[k]].insideComputeDomain) {
                dWR = ((FVM) elems[TT[k]]).volumeDerivative(null, u, elemData);
            } else {
                System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
            }
            if (elems[TT[k]].insideComputeDomain) {
                double[] WRp = elems[TT[k]].W;
                double sigmaR = ((FVM) elems[TT[k]]).FVMlimiter(dWR, FVMlimiter);
                for (int j = 0; j < nEqs; j++) {
                    double dW = 0;
                    for (int d = 0; d < dim; d++) {
                        dW = dW + (Int.faces[k].coordinatesFace[0][d] - elems[TT[k]].Xs[d]) * dWR[nEqs * d + j];
                    }
                    WR[j] = WRp[j] + sigmaR * dW;
                }
            } else {
                System.arraycopy(elems[TT[k]].W, 0, WR, 0, nEqs);
            }

        } else {
            WR = eqn.boundaryValue(WL, n[k][0], TT[k], elemData);
            System.arraycopy(dWL, 0, dWR, 0, dim * nEqs);
        }

        // inviscid flux in integration point
        double vn = 0;
        double[] Wale = new double[nEqs];
        if (eqn.isConvective()) {
            for (int d = 0; d < dim; d++) {
                vn = vn + u[d] * n[k][0][d];
            }
            if (TT[k] > -1) {
                for (int j = 0; j < nEqs; j++) {
                    Wale[j] = (WL[j] + WR[j]) / 2;
                }
            } else {
                System.arraycopy(WR, 0, Wale, 0, nEqs);
            }
            fn = eqn.numericalConvectiveFlux(WL, WR, n[k][0], TT[k], elemData);	// nevazky tok
        }

        // viscid flux in integration point
        if (eqn.isDiffusive()) {
            double[] Wc = new double[nEqs];
            double[] dWc = new double[nEqs * dim];
            for (int m = 0; m < nEqs; m++) {
                if (TT[k] > -1) {
                    Wc[m] = (WL[m] + WR[m]) / 2;
                } else {
                    Wc[m] = WR[m];
                }
                for (int d = 0; d < dim; d++) {
                    dWc[nEqs * d + m] = (dWL[nEqs * d + m] + dWR[nEqs * d + m]) / 2;
                }
            }
            fvn = eqn.numericalDiffusiveFlux(Wc, dWc, n[k][0], TT[k], elemData);
        }

        for (int m = 0; m < nEqs; m++) {
            if (eqn.isConvective()) {
                K[m] -= S[k] * (fn[m] - vn * Wale[m]);
            }
            if (eqn.isDiffusive()) {
                K[m] += S[k] * fvn[m];
            }
            if (KR != null) {
                if (eqn.isConvective()) {
                    KR[m] -= S[k] * (fn[m] - vn * Wale[m]);
                }
                if (eqn.isDiffusive()) {
                    KR[m] += S[k] * fvn[m];
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
            if (TT[k] > -1) {
                double[] WR = elems[TT[k]].W;
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

    public void computeMassMatrixAndBasisWeights() { // funkce pro generovani matic
        M[0][0] = area;
        Mo[0][0] = area;
        Mo2[0][0] = area;
        Is[0] = 1;
    }

    public void recalculateMassMatrix() {
        M[0][0] = area;
    }
    
    public boolean isJacobiMatrixAssembly(){
        return false;
    }
    
    public void limiter(boolean active){}
    
    public void residuumJacobi(double[][] ADiag, Implicit.Neighbour[] Sous){}
    
    public void residuumWallJacobi(int k, double[][] ADiag, Implicit.Neighbour Sous){}
}
