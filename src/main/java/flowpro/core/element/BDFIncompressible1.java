/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author obublik
 */
public class BDFIncompressible1 extends Implicit {
    double[] a1, a2;
    double beta;
    
    BDFIncompressible1(){
        super();
    }

    public int getOrder(){
        return 1;
    }
    
    public void init(FlowProProperties props){
        super.initImplicit();
        a1 = new double[nEqs];
        a2 = new double[nEqs];
        
        beta = 0;
        if (props.containsKey("beta")) {
            try {
                beta = props.getDouble("beta");
            } catch (IOException ex) {
                Logger.getLogger(BDFIncompressible1.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    void assembleRHS(double[] Rw, double[] a1, double[] a2) {
        double[][] M = elem.M;
        double[][] Mo = elem.Mo;
        double[] W = elem.W;
        double[] Wo = elem.Wo;

        System.arraycopy(Rw, 0, RHS_loc, 0, Rw.length);
        for (int m = 0; m < nEqs; m++) {
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    RHS_loc[nBasis * m + i] -= M[i][j] * a1[m] * W[m * nBasis + j] + Mo[i][j] * a2[m] * Wo[m * nBasis + j];
                }
            }
        }
    }

    // Generovani radku globalni matice a vektoru prave strany
    public void assembleJacobiMatrix(double dt, double dto) {
        // coefs
        a1[0] = beta / dt;
        a2[0] = -beta / dt;
        for(int m = 1; m < nEqs; m++){
            a1[m] = 1 / dt;
            a2[m] = -1 / dt;
        }
        
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
        assembleRHS(Rw, a1, a2);

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
                    ADiag[nBasis * m + i][nBasis * m + j] += a1[m] * M[i][j];
                }
            }
        }
    }
}
