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

    int[] Iblock;
    int[] Jblock;
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
        int dofs = A.getDofs();
        int nnz = A.getNNZ(); 
        int[] Icoo = A.getRowIndexesCOO();
        int[] Jcoo = A.getColumnIndexesCOO();
        Iin = new boolean[dofs];
        Jin = new boolean[dofs];
        for(int i = 0; i < nIBlock; i++){
            Iin[Iblock[i]] = true;
        }
        for(int i = 0; i < nJBlock; i++){
            Jin[Jblock[i]] = true;
        }
        int s = 0;
        for (int i = 0; i < nnz; i++) { // get block nnz
            if(Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                s++;
            }
        }
        IblockCoo = new int[s];
        JblockCoo = new int[s];
        Hblock = new double[s];
        s = 0;
        for (int i = 0; i < nnz; i++) { // fill block indexes
            if(Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                IblockCoo[s] = Icoo[i];
                JblockCoo[s] = Jcoo[i];
                s++;
            }
        }
    }
    
    void setH(){
        double[] H = A.getData();
        int s = 0;
        int[] Icoo = A.getRowIndexesCOO();
        int[] Jcoo = A.getColumnIndexesCOO();
        for (int i = 0; i < A.getNNZ(); i++) { // fill block indexes
            if(Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                Hblock[s] = H[i];
                s++;
            }
        }
    }
}
