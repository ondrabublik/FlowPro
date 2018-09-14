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
    double[] Wo;
    double[][] MWHistory;
    double[][] RHSHistory;
    double[] timeHistory;
    int nEqs;
    int nBasis;

    StorageTimeLevels(int timeLevels, int nBasis, int nEqs) {
        this.timeLevels = timeLevels;
        this.nBasis = nBasis;
        this.nEqs = nEqs;

    }

    double[] init(double t, double[] initW, double[][] M, Basis basis) {
        this.timeLevels = timeLevels;
        size = nBasis * nEqs;
        double[] W = new double[size];
        Wo = new double[size];
        MWHistory = new double[timeLevels][size];
        RHSHistory = new double[timeLevels][size];
        timeHistory = new double[timeLevels];
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
            for (int k = 0; k < timeLevels; k++) {
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

        for (int k = 0; k < timeLevels; k++) {
            timeHistory[k] = t;
        }

        position = timeLevels - 1;

        return W;
    }

    void calculateRHS(double[] RHS, double[][] M, double[] W, double[] WCoef, double[] RHSCoef) {
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    RHS[nBasis * m + i] -= WCoef[0] * M[i][j] * W[nBasis * m + j];
                }
                for (int k = 1; k < WCoef.length; k++) {
                    RHS[nBasis * m + i] -= WCoef[k] * MWHistory[(position - k) % timeLevels][nBasis * m + i];
                }
                for (int k = 0; k < RHSCoef.length; k++) {
                    RHS[nBasis * m + i] += RHSCoef[k] * RHSHistory[(position - k)%timeLevels][nBasis * m + i];
                } 
           }
        }
    }

    // ulozeni W do Wo
    void saveTimeLevel(double t, double[][] M, double[] W, double[] RHS) {
        position = plus(position);
        if (RHS != null) { //implicit method
            for (int m = 0; m < nEqs; m++) {
                for (int i = 0; i < nBasis; i++) {
                    MWHistory[position][nBasis * m + i] = 0;
                    for (int j = 0; j < nBasis; j++) {
                        MWHistory[position][nBasis * m + i] += M[i][j] * W[nBasis * m + j];
                    }
                }
            }
            timeHistory[position] = t;
            System.arraycopy(RHS, 0, RHSHistory[position], 0, size);
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

    double getDto(double t){
        return t - timeHistory[position];
    }
    
    int plus(int i) {
        return (i++) % timeLevels;
    }

    int minus(int i) {
        return (i--) % timeLevels;
    }
}
