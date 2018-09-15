/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core;

import flowpro.core.basis.Basis;

/**
 *
 * @author obublik
 */
public class StorageTimeLevels {

    int size;
    int posW;
    int posRHS;
    double[] Wo;
    double[][] MWHistory;
    double[][] RHSHistory;
    int nEqs;
    int nBasis;
    int nMW;
    int nRHS;

    StorageTimeLevels(int nW, int nRHS, int nBasis, int nEqs) {
        this.nBasis = nBasis;
        this.nEqs = nEqs;
        this.nMW = nW;
        this.nRHS = nRHS;
    }

    double[] init(double[] initW, double[][] M, Basis basis) {
        size = nBasis * nEqs;
        double[] W = new double[size];
        Wo = new double[size];
        if (nMW > 0) {
            MWHistory = new double[nMW][size];
        }
        if (nRHS > 0) {
            RHSHistory = new double[nRHS][size];
        }
        this.nEqs = nEqs;
        this.nBasis = nBasis;

        // fill the solution vector with initial condition
        if (initW.length == size) {
            System.arraycopy(initW, 0, W, 0, size);
            System.arraycopy(initW, 0, Wo, 0, size);
        } else {
            switch (basis.basisType) {
                case "lagrange":
                    for (int i = 0; i < nBasis; i++) {
                        for (int j = 0; j < nEqs; j++) {
                            W[j * nBasis + i] = initW[j];
                            Wo[j * nBasis + i] = initW[j];
                        }
                    }
                    break;
                case "orthogonal":
                case "taylor":
                    for (int j = 0; j < nEqs; j++) {
                        W[j * nBasis] = initW[j];
                        Wo[j * nBasis] = initW[j];
                    }
                    break;
            }
        }

        // init MW history
        if (M != null) {
            for (int k = 0; k < nMW; k++) {
                for (int m = 0; m < nEqs; m++) {
                    for (int i = 0; i < nBasis; i++) {
                        MWHistory[k][nBasis * m + i] = 0;
                        for (int j = 0; j < nBasis; j++) {
                            MWHistory[k][nBasis * m + i] += M[i][j] * W[nBasis * m + j];
                        }
                    }
                }
            }
        }

        posW = nMW - 1;

        return W;
    }

    void calculateRHS(double[] RHS, double[][] M, double[] W, double[] WCoef, double[] RHSCoef) {
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    RHS[nBasis * m + i] -= WCoef[0] * M[i][j] * W[nBasis * m + j];
                }
                if (nMW > 0) {
                    for (int k = 1; k < WCoef.length; k++) {
                        RHS[nBasis * m + i] -= WCoef[k] * MWHistory[(posW - k + 1 + nMW) % nMW][nBasis * m + i];
                    }
                }
                if (nRHS > 0) {
                    for (int k = 1; k < RHSCoef.length; k++) {
                        RHS[nBasis * m + i] += RHSCoef[k] * RHSHistory[(posRHS - k + 1 + nRHS) % nRHS][nBasis * m + i];
                    }
                }
            }
        }
    }

    // ulozeni W do Wo
    void storeTimeLevel(double[][] M, double[] W, double[] RHS) {
        posW = (posW+1) % nMW;
        posRHS = (posRHS++) % nRHS;
        if (RHS != null) { //implicit method
            if (nMW > 0) {
                for (int m = 0; m < nEqs; m++) {
                    for (int i = 0; i < nBasis; i++) {
                        MWHistory[posW][nBasis * m + i] = 0;
                        for (int j = 0; j < nBasis; j++) {
                            MWHistory[posW][nBasis * m + i] += M[i][j] * W[nBasis * m + j];
                        }
                    }
                }
            }
            if (nRHS > 0) {
                System.arraycopy(RHS, 0, RHSHistory[posRHS], 0, size);
            }
        }
        System.arraycopy(W, 0, Wo, 0, size);
    }

    double calculateResiduumW(double dt, double[] W) {
        double rez = 0;
        for (int i = 0; i < size; i++) {
            if (dt > 0) {
                rez = rez + Math.abs(W[i] - Wo[i]) / dt;
            }
        }
        return rez;
    }

    // ulozeni Wo do W (pri potizich s resenim)
    void getTimeLevelBack(double[] W) {
        System.arraycopy(Wo, 0, W, 0, size);
    }
}
