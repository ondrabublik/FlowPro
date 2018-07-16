package flowpro.core.LinearSolvers;

import eigenwrapper.EigenWrapper;
import flowpro.core.Mesh;
import flowpro.core.Parameters;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umfpackwrapper.UmfpackWrapper;

/**
 *
 * @author obublik
 */
public class ExternCppSolver extends LinearSolver {

    int dofs;
    int nnz;
    public int[] IA, JA;
    public double[] HA, b; // CRS matrix format
    //EigenWrapper eigenWrap;
    UmfpackWrapper umfWrap;

    ExternCppSolver(Mesh.Element[] elems, int dofs, Parameters par) {
        this.elems = elems;
        this.dofs = dofs;
        nnz = getNNZ();
        IA = new int[nnz];
        JA = new int[nnz];
        HA = new double[nnz];
        b = new double[dofs];
        buildCRSStructure();
        //eigenWrap = new EigenWrapper();
        umfWrap = new UmfpackWrapper();
    }

    @Override
    public boolean solve(double[] x) {
        try {
            fillAb();
            umfWrap.solveUmfpack(IA, JA, HA, b, x);
            //buildMatrix();
            //eigenWrap.solveEigen(IA, JA, HA, b, x);
            return true;
        } catch (Exception ex) {
            Logger.getLogger(ExternSolver.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }
    }

    public void buildMatrix() throws IOException {
        int s = 0;
        for (Mesh.Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            int[] glob = elem.gi_U;
            double[][] Ad = elem.ADiag;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    IA[s] = glob[i];
                    JA[s] = glob[j];
                    HA[s] = Ad[j][i];
                    s++;
                }
            }
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    int[] globe = elems[elem.TT[k]].gi_U;
                    double[][] An = elem.ANeighs[k].MR;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < ne; j++) {
                            IA[s] = glob[i];
                            JA[s] = globe[j];
                            HA[s] = An[j][i];
                            s++;
                        }
                    }
                }
            }
        }
        
        s = 0;
        for (Mesh.Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            for (int i = 0; i < n; i++) {
                b[s] = elem.RHS_loc[i];
                s++;
            }
        }
    }
    
    public void buildCRSStructure() {
        int[] IAaux = new int[nnz];
        int si = 0;
        int sj = 0;
        for (Mesh.Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            int[] glob = elem.gi_U;
            for (int i = 0; i < n; i++) { // row cycle
                for (int j = 0; j < n; j++) { // column cycle
                    JA[sj] = glob[j];
                    sj++;
                }
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TT[k] > -1) {
                        int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                        int[] globe = elems[elem.TT[k]].gi_U;
                        for (int j = 0; j < ne; j++) {
                            JA[sj] = globe[j];
                            sj++;
                        }
                    }
                }
                IAaux[si] = sj;
                si++;
            }
        }
        IA = new int[si];
        System.arraycopy(IAaux, 0, IA, 0, si);
        IAaux = null;
    }
    
    public void fillAb() {
        int sj = 0;
        for (Mesh.Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            double[][] Ad = elem.ADiag;
            for (int i = 0; i < n; i++) { // row cycle
                for (int j = 0; j < n; j++) { // column cycle
                    HA[sj] = Ad[j][i];
                    sj++;
                }
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TT[k] > -1) {
                        int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                        double[][] An = elem.ANeighs[k].MR;
                        for (int j = 0; j < ne; j++) {
                            HA[sj] = An[j][i];
                            sj++;
                        }
                    }
                }
            }
        }
        
        int s = 0;
        for (Mesh.Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            for (int i = 0; i < n; i++) {
                b[s] = elem.RHS_loc[i];
                s++;
            }
        }
    }

    public int getNNZ() {
        int s = 0;
        for (Mesh.Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            s += n * n;
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    s += n * ne;
                }
            }
        }
        return s;
    }
}