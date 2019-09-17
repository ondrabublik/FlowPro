package flowpro.core.element;

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
