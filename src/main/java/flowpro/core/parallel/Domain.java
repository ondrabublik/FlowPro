package flowpro.core.parallel;

import flowpro.core.curvedBoundary.FaceCurvature;
import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class Domain implements Serializable {

    public final int nDoms;
    public final int nElems;
    public final int nPoints;

    private final Subdomain[] subdoms;

    public Domain(int[] elemsOrder, int[] elemsType, int[][] TT, int[][] TP, int[][] TEale, int[][] TEshift, FaceCurvature[] fCurv, double[][] initW, int[] metisMap, int nDoms, int nOverlap, int nPoints) throws IOException {
        this.nDoms = nDoms;

        /*
         if s = dAux[i][j] takes the value
         - s = 0 then the i-th element is not in the j-th domain,
         - s = 1 then the i-th element is in the interior of the j-th domain (not in the overlap),
         - s > 1 then the i-th element is in the (s-1)-st layer of the overlap.
         */
        nElems = TT.length;
        this.nPoints = nPoints;
        
        int[][] dAux = new int[nElems][nDoms];
        // identify interior elements
        for (int i = 0; i < nElems; i++) {
            dAux[i][metisMap[i]] = 1;
        }

        // identify overlaping elements
        for (int j = 0; j < nDoms; j++) {
            for (int p = 1; p < nOverlap + 1; p++) {
                for (int i = 0; i < nElems; i++) {
                    if (dAux[i][j] == p) {
                        int nEdges = TT[i].length;
                        for (int k = 0; k < nEdges; k++) {
                            if (TT[i][k] > -1 && dAux[TT[i][k]][j] == 0) {
                                dAux[TT[i][k]][j] = p + 1;
                            }
                        }
                    }
                }
            }
        }

        // map generation
        subdoms = new Subdomain[nDoms];
        for (int j = 0; j < nDoms; j++) {
            int nSub = 0; // pocet prvku v podoblasti
            for (int i = 0; i < nElems; i++) {
                if (dAux[i][j] != 0) {
                    nSub++;
                }
            }
            int[] mapL2G = new int[nSub]; // local to global indexing
            int[] mapG2L = new int[nElems]; // global to local indexing
            Arrays.fill(mapG2L, -100);
            int s = 0;
            for (int i = 0; i < nElems; i++) {
                if (dAux[i][j] != 0) {
                    mapL2G[s] = i;
                    mapG2L[i] = s;
                    s++;
                }
            }

            // local mesh and indexing creation
            int[][] subTP = new int[nSub][];
            for (int i = 0; i < nSub; i++) {
                int nPoint = TP[mapL2G[i]].length;
                subTP[i] = new int[nPoint];              
                System.arraycopy(TP[mapL2G[i]], 0, subTP[i], 0, nPoint);
            }
            
            // local mesh and indexing creation
            int[] subElemsOrder = new int[nSub];
            int[] subElemsType = new int[nSub];
            int[][] subTT = new int[nSub][];
            int[][] subTEale = new int[nSub][];
            int[][] subTEshift = new int[nSub][];
            FaceCurvature[] subfCurv = new FaceCurvature[nSub];
            double[][] subInitW = new double[nSub][];
            for (int i = 0; i < nSub; i++) {
                int nEdges = TT[mapL2G[i]].length;
                subElemsOrder[i] = elemsOrder[mapL2G[i]];
                subElemsType[i] = elemsType[mapL2G[i]];
                subfCurv[i] = fCurv[mapL2G[i]];
                subTT[i] = new int[nEdges];
                subTEale[i] = new int[nEdges];
                subTEshift[i] = new int[nEdges];
                for (int k = 0; k < nEdges; k++) {
                    subTT[i][k] = TT[mapL2G[i]][k];
                    subTEale[i][k] = TEale[mapL2G[i]][k];
                    subTEshift[i][k] = TEshift[mapL2G[i]][k];
                    if (subTT[i][k] > -1) {
                        subTT[i][k] = mapG2L[TT[mapL2G[i]][k]];
                    }
                }
                int nEqBas = initW[mapL2G[i]].length;
                subInitW[i] = new double[nEqBas];
                System.arraycopy(initW[mapL2G[i]], 0, subInitW[i], 0, nEqBas);
            }

            // count number of boudary, load and save elements
            int nBoundaryElems = 0;
            int nLoadElems = 0;
            int nSaveElems = 0;
            int nInterior = 0;
            for (int loc = 0; loc < nSub; loc++) {
                int glob = mapL2G[loc];
                if (dAux[glob][j] == nOverlap + 1) {
                    nBoundaryElems++;
                }
                if (dAux[glob][j] > 1) {
                    nLoadElems++;
                }
                if (dAux[glob][j] == 1) {
                    nInterior++;
                }
                if (dAux[glob][j] == 1) {
                    for (int p = 0; p < nDoms; p++) {
                        if (dAux[glob][p] > 1) {
                            nSaveElems++;
                        }
                    }
                }
            }

            int[] boundary = new int[nBoundaryElems];
            int[] load = new int[nLoadElems];
            int[] save = new int[nSaveElems];
            int[] interior = new int[nInterior];

            int boundaryIdx = 0;
            int loadIdx = 0;
            int saveIdx = 0;
            int interiorIdx = 0;
            for (int loc = 0; loc < nSub; loc++) {
                int glob = mapL2G[loc];

                // local indexes of the boundary elements
                if (dAux[glob][j] == nOverlap + 1) {
                    boundary[boundaryIdx] = loc;
                    boundaryIdx++;
                }

                // global indexes of elements which are to be loaded from the central structure
                if (dAux[glob][j] > 1) {
                    load[loadIdx] = glob;
                    loadIdx++;
                }

                if (dAux[glob][j] == 1) {
                    interior[interiorIdx] = loc;
                    interiorIdx++;
                }

                // local indexes of elements which are to be saved into the central structure
                if (dAux[glob][j] == 1) {
                    for (int p = 0; p < nDoms; p++) {
                        if (dAux[glob][p] > 1) {
                            save[saveIdx] = loc;
                            saveIdx++;
                            break;
                        }
                    }
                }
            }

            subdoms[j] = new Subdomain(nSub, mapL2G, mapG2L, subElemsOrder, subElemsType, subTP, subTT, subTEale, subTEshift, subfCurv, subInitW, boundary, load, save, interior);
        }
    }

    public Subdomain getSubdomain(int n) {
        return subdoms[n];
    }

    public class Subdomain implements Serializable {

        public final int nElems;
        public final int[] mapL2G;
        public final int[] mapG2L;
        public final int[] elemsOrder;
        public final int[] elemsType;
        public final int[][] TP;
        public final int[][] TT;
        public final int[][] TEale;
        public final int[][] TEshift;
        public final FaceCurvature[] fCurv;
        public final double[][] initW;
        public final int[] boundary;
        public final int[] load;
        public final int[] save;
        public final int[] interior;
        
        /**
         *
         * @param nElems
         * @param mapL2G local to global indexes map
         * @param mapG2L global to local indexes map
         * @param TP
         * @param TT
         * @param TEale
         * @param boundary indexes of boundary elements
         * @param load global indexes of elements located in the outer overlap
         * @param save local indexes of elements located in the inner overlap
         * @param interior local indexes of elements in the interior (not in the
         * outer overlap)
         */
        private Subdomain(int nElems, int[] mapL2G, int[] mapG2L, int[] elemsOrder, int[] elemsType, int[][] TP, int[][] TT, int[][] TEale, int[][] TEshift, FaceCurvature[] fCurv, double[][] initW, int[] boundary,
                int[] load, int[] save, int[] interior) {
            this.nElems = nElems;
            this.mapL2G = mapL2G;
            this.mapG2L = mapG2L;
            this.elemsOrder = elemsOrder;
            this.elemsType = elemsType;
            this.TP = TP;
            this.TT = TT;
            this.TEale = TEale;
            this.TEshift = TEshift;
            this.fCurv = fCurv;
            this.initW = initW;
            this.boundary = boundary;
            this.load = load;
            this.save = save;
            this.interior = interior;
        }
    }
}
