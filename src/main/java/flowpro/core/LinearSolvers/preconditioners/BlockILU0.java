/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.api.Mat;
import flowpro.core.LinearSolvers.SparseMatrix;
import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;

/**
 *
 * @author obublik
 */
public class BlockILU0 extends Preconditioner {

    SparseMatrix A;
    int n, nThreads;
    int[][] indexes;
    double[][][] diagonalInverse;
    double[][][][] ILU;
    Element[] elems;
    int nElem;

    BlockILU0(Parameters par) {
        nThreads = par.nThreads;
    }

    @Override
    public void setMatrix(SparseMatrix A) {
        elems = A.getElems();
        nElem = elems.length;
        diagonalInverse = new double[nElem][][];
        ILU = new double[nElem][][][];
        indexes = new int[nElem][];
        for(int i = 0; i < nElem; i++){
            int nNeigh = 0;
            for(int k = 0; k < elems[i].nFaces; k++){
                if(elems[k].TT[k] > -1){
                    nNeigh++;
                }
            }
            ILU[i] = new double[nNeigh][][];
            indexes[i] = new int[nNeigh];
            int s = 0;
            int ne = elems[i].nBasis*elems[i].getNEqs();
            for(int k = 0; k < elems[i].nFaces; k++){
                if(elems[k].TT[k] > -1){
                    int me = elems[elems[i].TT[k]].nBasis*elems[elems[i].TT[k]].getNEqs();
                    ILU[i][s] = new double[ne][me];
                    //indexes[i][s] = 
                }
            }
        }
    }

    @Override
    public void factor() {
        for(int i = 0; i < nElem; i++){
            diagonalInverse[i] = Mat.invert(elems[i].ADiag);
        }
//        for (int i = 1; i < n; i++) {
//            for (int k = 0; k < i - 1; k++) {
//                A(i, k) = A(i, k) / A(k, k);
//                for (int j = k + 1; j < n; j++) {
//                    A(i, j) = A(i, j) - A(i, k) * A(k, j);
//                }
//            }
//        }
    }

    @Override
    public void apply(double[] x, double[] b) {

    }
}
