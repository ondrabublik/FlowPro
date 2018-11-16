/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

/**
 *
 * @author obublik
 */
public class SparseMatrixBlock {

    int nnz, dofs;
    int[] Iblock, Imap;
    int[] Jblock, Jmap;
    boolean[] Iin, Jin;
    int nIBlock;
    int nJBlock;
    SparseMatrix A;
    int[] IblockCoo;
    int[] JblockCoo;
    double[] Hblock;

    SparseMatrixBlock(int[] Iblock, int[] Jblock, SparseMatrix A) {
        this.Iblock = Iblock;
        this.Jblock = Jblock;
        this.A = A;
        nIBlock = Iblock.length;
        nJBlock = Jblock.length;

        extractBlock();
    }

    final void extractBlock() {
        int dofsA = A.getDofs();
        int nnzA = A.getNNZ();
        int[] Icoo = A.getRowIndexesCOO();
        int[] Jcoo = A.getColumnIndexesCOO();
        Iin = new boolean[dofsA];
        Jin = new boolean[dofsA];
        Imap = new int[dofsA];
        Jmap = new int[dofsA];
        for (int i = 0; i < nIBlock; i++) {
            Iin[Iblock[i]] = true;
            Imap[Iblock[i]] = i;
        }
        for (int i = 0; i < nJBlock; i++) {
            Jin[Jblock[i]] = true;
            Jmap[Jblock[i]] = i;
        }
        int s = 0;
        for (int i = 0; i < nnzA; i++) { // get block nnzA
            if (Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                s++;
            }
        }
        nnz = s;
        IblockCoo = new int[nnz];
        JblockCoo = new int[nnz];
        Hblock = new double[nnz];
        
        double[] H = A.getData();
        s = 0;
        for (int i = 0; i < nnzA; i++) { // fill block indexes
            if (Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                IblockCoo[s] = Imap[Icoo[i]];
                JblockCoo[s] = Jmap[Jcoo[i]];
                s++;
            }
        }
    }

    void setH() {
        double[] H = A.getData();
        int s = 0;
        int[] Icoo = A.getRowIndexesCOO();
        int[] Jcoo = A.getColumnIndexesCOO();
        for (int i = 0; i < A.getNNZ(); i++) { // fill block indexes
            if (Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                Hblock[s] = H[i];
                s++;
            }
        }
    }
    
    int getNRows() {
        return nIBlock;
    }

    int getNColumn() {
        return nJBlock;
    }
    
    int getNNz() {
        return nnz;
    }

    public void SubstrMult(double[] y, double[] b, double[] x) {
        for(int i = 0; i < nIBlock; i++){
            y[Iblock[i]] = b[Iblock[i]];
        }
        for (int i = 0; i < IblockCoo.length; i++) {
            y[IblockCoo[i]] -= Hblock[i] * x[JblockCoo[i]];
        }
    }
}
