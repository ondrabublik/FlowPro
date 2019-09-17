/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.Mat;

/**
 *
 * @author obublik
 */
public class ImplicitTimeIntegration extends TimeIntegration {

    public double[][] ADiag;
    double[][] PrecondJacobi; // inverze ADiag
    public double[] RHS_loc;
    public Neighbour[] ANeighs; //matice sousedu
    public int[] TT;
    public Element[] elems;

    ImplicitTimeIntegration(Element elem) {
        super(elem);
        isImplicit = true;

        TT = elem.TT;
        ADiag = new double[nEqs * nBasis][nEqs * nBasis];
        RHS_loc = new double[nEqs * nBasis];
        elems = elem.elems;
    }

    public void init(){
        alocateNeigbourhsLinearSolver();
    }
    
    public void alocateNeigbourhsLinearSolver() {
        ANeighs = new Neighbour[nFaces];
        for (int k = 0; k < nFaces; k++) {
            if (TT[k] > -1) {
                ANeighs[k] = new Neighbour(TT[k], nBasis, elems[TT[k]].nBasis, nEqs);
            } else {
                ANeighs[k] = new Neighbour(TT[k], nBasis, 0, nEqs);
            }
        }
    }

    void assembleRHS(double[] Rw, double[] a1, double[] a2, double[] a3) {
        double[][] M = elem.M;
        double[][] Mo = elem.Mo;
        double[][] Mo2 = elem.Mo2;
        double[] W = elem.W;
        double[] Wo = elem.Wo;
        double[] Wo2 = elem.Wo2;

        System.arraycopy(Rw, 0, RHS_loc, 0, Rw.length);
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    RHS_loc[nBasis * m + i] -= M[i][j] * a1[m] * W[m * nBasis + j] + Mo[i][j] * a2[m] * Wo[m * nBasis + j] + Mo2[i][j] * a3[m] * Wo2[m * nBasis + j];
                }
            }
        }
    }

    // Generovani radku globalni matice a vektoru prave strany
    public void assembleJacobiMatrix(double[] a1, double[] a2, double[] a3, double[] dual) {

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
        elem.residuum(V, Rw, RwNeigh);
        // assemble rhs
        assembleRHS(Rw, a1, a2, a3);

//        if (elem.par.useJacobiMatrix && elem.eqn.isEquationsJacobian()) { // fast assemble when jacobian of equations is known
//            residuumWithJacobian(ADiag, ANeighs);
//        } else { // slow assemble when jacobian of equations is unknown
        double h = elem.par.h;
        for (int i = 0; i < nBasis * nEqs; i++) {
            for (int j = 0; j < Rw.length; j++) {
                ADiag[i][j] = -Rw[j];
            }
            V[i] = h;
            elem.residuum(V, ADiag[i], RwNeighH);
            V[i] = 0;
            for (int k = 0; k < nFaces; k++) {
                if (TT[k] > -1 && elems[TT[k]].insideComputeDomain) {
                    double[][] Aaux = ((ImplicitTimeIntegration) elems[TT[k]].ti).ANeighs[elem.faceIndexReverse[k]].A;
                    for (int j = 0; j < RwNeighH[k].length; j++) {
                        Aaux[i][j] = (RwNeighH[k][j] - RwNeigh[k][j]) / h;
                    }
                }
            }
        }
        Mat.divide(ADiag, -h);

        // pricteni matice hmotnosti
        double[][] M = elem.M;
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    ADiag[nBasis * m + i][nBasis * m + j] += (a1[m] + dual[m]) * M[i][j];
                }
            }
        }
    }


    public void updateW(double[] x) {
        int[] gi_U = elem.gi_U;
        double[] W = elem.W;
        for (int i = 0; i < nEqs * nBasis; i++) {
            W[i] += x[gi_U[i]];
        }
    }

    public void updateRHS(double[] x) {
        int[] gi_U = elem.gi_U;
        int n = nEqs * nBasis;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                RHS_loc[i] = RHS_loc[i] - ADiag[j][i] * x[gi_U[j]];
            }
            for (int k = 0; k < nFaces; k++) {
                if (TT[k] > -1) {
                    for (int j = 0; j < nEqs * elems[TT[k]].nBasis; j++) {
                        RHS_loc[i] = RHS_loc[i] - ANeighs[k].A[j][i] * x[elems[TT[k]].gi_U[j]];
                    }
                }
            }
        }
    }

    public double sqr() {
        double nrm = 0;
        for (int i = 0; i < nEqs * nBasis; i++) {
            nrm = nrm + RHS_loc[i] * RHS_loc[i];
        }
        return nrm;
    }
}
