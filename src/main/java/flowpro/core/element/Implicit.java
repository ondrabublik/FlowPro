/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

/**
 *
 * @author obublik
 */
public abstract class Implicit extends TimeIntegrationElement {
    
    public double[][] ADiag;
    public Neighbour[] ANeighs; //matice sousedu
    public double[] RHS_loc;
    public int[] TT;
    public Element[] elems;
    
    Implicit() {
        super();
    }
    
    abstract public void assembleJacobiMatrix(double dt, double dto);
    
    public boolean isImplicit(){
        return true;
    }
    
    public String getLocalSolverType(){
        return "localimplicit";
    }
    
    public void init(){
        isImplicit = true;
        TT = elem.TT;
        ADiag = new double[nEqs * nBasis][nEqs * nBasis];
        RHS_loc = new double[nEqs * nBasis];
        elems = elem.elems;
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
    
    public class Neighbour {

    public double[][] A;
    double[] V;  // e_j * h
    int nr, nBasis, neR, typ; // neR - pocet bazovych funkci souseda

    public Neighbour(int typ, int nBasis, int neR, int nr) {
        if (typ > -1) {
            A = new double[nr * neR][nr * nBasis];
        }
        V = new double[neR * nr];
        this.nr = nr;
        this.nBasis = nBasis;
        this.neR = neR;
        this.typ = typ;
    }

    void nullMatrix() {
        if (typ > -1) {
            for (int i = 0; i < nBasis * nr; i++) {
                for (int j = 0; j < neR * nr; j++) {
                    A[j][i] = 0;
                }
            }
        }
    }
}
}
