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
        int dofs = A.getDofs();
        int nnz = A.getNNZ();
        int[] Icoo = A.getRowIndexesCOO();
        int[] Jcoo = A.getColumnIndexesCOO();
        Iin = new boolean[dofs];
        Jin = new boolean[dofs];
        Imap = new int[dofs];
        Jmap = new int[dofs];
        for (int i = 0; i < nIBlock; i++) {
            Iin[Iblock[i]] = true;
            Imap[Iblock[i]] = i;
        }
        for (int i = 0; i < nJBlock; i++) {
            Jin[Jblock[i]] = true;
            Jmap[Jblock[i]] = i;
        }
        int s = 0;
        for (int i = 0; i < nnz; i++) { // get block nnz
            if (Iin[Icoo[i]] && Jin[Jcoo[i]]) {
                s++;
            }
        }
        IblockCoo = new int[s];
        JblockCoo = new int[s];
        Hblock = new double[s];
        s = 0;
        for (int i = 0; i < nnz; i++) { // fill block indexes
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

    public void SubstrMult(double[] y, double[] b, double[] x) {
        System.arraycopy(b, 0, y, 0, nIBlock);
        for (int i = 0; i < IblockCoo.length; i++) {
            y[i] -= Hblock[i] * x[JblockCoo[i]];
        }
    }
}
