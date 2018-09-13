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

    int timeLevels;
    int size;
    int position;
    double[][] WHistory;
    double[][] MWHistory;
    double[][] RHSHistory;
    double[] timeHistory;
    int nEqs;
    int nBasis;

    StorageTimeLevels(int timeLevels, double[] initW, int nBasis, int nEqs, Basis basis) {
        this.timeLevels = timeLevels;
        size = nBasis * nEqs;
        WHistory = new double[timeLevels][size];
        MWHistory = new double[timeLevels][size];
        RHSHistory = new double[timeLevels][size];
        timeHistory = new double[timeLevels];
        this.nEqs = nEqs;
        this.nBasis = nBasis;

        // fill the solution vector with initial condition
        if (initW.length == size) {
            for (int k = 0; k < timeLevels; k++) {
                System.arraycopy(initW, 0, WHistory[k], 0, size);
            }
        } else {
            switch (basis.basisType) {
                case "lagrange":
                    for (int k = 0; k < timeLevels; k++) {
                        for (int i = 0; i < nBasis; i++) {
                            for (int j = 0; j < nEqs; j++) {
                                WHistory[k][j * nBasis + i] = initW[j];
                            }
                        }
                    }
                    break;
                case "orthogonal":
                case "taylor":
                    for (int k = 0; k < timeLevels; k++) {
                        for (int j = 0; j < nEqs; j++) {
                            WHistory[k][j * nBasis] = initW[j];
                        }
                    }
                    break;
            }
        }
        position = timeLevels - 1;
    }

    void calculateRHS(double[] RHS, double[][] M, double[] W, double[] coef) {
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    RHS[nBasis * m + i] -= M[i][j] * coef[0] * W[m * nBasis + j];
                }
                for (int k = 1; k < coef.length; k++) {
                    RHS[nBasis * m + i] -= coef[k] * MWHistory[(position + k)%timeLevels][nBasis * m + i];
                }
//                for (int k = 0; k < coef.length; k++) {
//                    RHS[nBasis * m + i] += coef[k] * RHSHistory[(position + k)%timeLevels][nBasis * m + i];
//                }
            }
        }
    }

    // ulozeni W do Wo
    void saveTimeLevel(double t, double[] W, double[] MW, double[] RHS) {
        position = plus(position);
        timeHistory[position] = t;
        System.arraycopy(W, 0, WHistory[position], 0, size);
        System.arraycopy(MW, 0, MWHistory[position], 0, size);
        System.arraycopy(RHS, 0, RHSHistory[position], 0, size);
    }

    double calculateResiduumW(double t, double[] W) {
        double rez = 0;
        double dt = t - timeHistory[position];
        for (int i = 0; i < size; i++) {
            if (dt > 0) {
                rez = rez + Math.abs(W[i] - WHistory[minus(position)][i]) / dt;
            }
        }
        return rez;
    }

    // ulozeni Wo do W (pri potizich s resenim)
    void getTimeLevelBack(double[] W) {
        System.arraycopy(WHistory[position], 0, W, 0, size);
    }

    int plus(int i) {
        return (i++) % timeLevels;
    }

    int minus(int i) {
        return (i--) % timeLevels;
    }
}
