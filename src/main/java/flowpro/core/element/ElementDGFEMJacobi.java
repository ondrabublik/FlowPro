///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package flowpro.core.element;
//
//import flowpro.api.Mat;
//import flowpro.core.Mesh;
//import flowpro.core.curvedBoundary.FaceCurvature;
//import flowpro.core.elementType.ElementType;
//import java.io.IOException;
//
///**
// *
// * @author obublik
// */
//public class ElementDGFEMJacobi extends Element {
//    
//    public ElementDGFEMJacobi(int index, double[][] vertices, double[][] Uinit, double[] wallDistance, double[][] externalField, int[] TT, int[] TP, int[] TEale, int[] TEshift, double[][] shift, FaceCurvature fCurv, double[][] blendFun, double[] initW,
//            Mesh mesh, ElementType elemType) throws IOException {
//        super(index, vertices, Uinit, wallDistance, externalField, TT, TP, TEale, TEshift, shift, fCurv, blendFun, initW, mesh, elemType);
//        
//    }
//    
//    public void residuum(double[][] ADiag, Neighbour[] Sous) {
//        // vypocet toku hranici
//        for (int k = 0; k < nFaces; k++) { // opakovani pres jednotlive steny
//            residuumWall(k, ADiag, Sous[k]);
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
//    public void residuumWall(int k, double[][] ADiag, Neighbour Sous) {
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
//}