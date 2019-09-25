/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.FlowProProperties;
import flowpro.api.Mat;

/**
 *
 * @author obublik
 */
public class BDF2 extends Implicit {
    
    boolean useJacobiMatrixForAssembly;

    BDF2() {
        super();
    }

    public int getOrder() {
        return 2;
    }

    public void init(FlowProProperties props) {
        super.initImplicit();
        useJacobiMatrixForAssembly = elem.isJacobiMatrixAssembly();
    }

    void assembleRHS(double[] Rw, double a1, double a2, double a3) {
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
                    RHS_loc[nBasis * m + i] -= M[i][j] * a1 * W[m * nBasis + j] + Mo[i][j] * a2 * Wo[m * nBasis + j] + Mo2[i][j] * a3 * Wo2[m * nBasis + j];
                }
            }
        }
    }

    // Generovani radku globalni matice a vektoru prave strany
    public void assembleJacobiMatrix(double dt, double dto) {
        // coefs
        double a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt);
        double a2 = -(dt + dto) / (dt * dto);  // -2/dt;
        double a3 = dt / (dto * (dt + dto));  // 1/(2*dt);

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

        if (useJacobiMatrixForAssembly) { // fast assemble when jacobian of equations is known
            ADiag = new double[nEqs * nBasis][nEqs * nBasis];
            for (Neighbour neigh : ANeighs) {
                neigh.nullMatrix();
            }
            elem.residuumJacobi(ADiag, ANeighs);
        } else { // slow assemble when jacobian of equations is unknown
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
                        double[][] Aaux = ((Implicit) elems[TT[k]].ti).ANeighs[elem.faceIndexReverse[k]].A;
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
                        ADiag[nBasis * m + i][nBasis * m + j] += a1 * M[i][j];
                    }
                }
            }
        }
    }
}
