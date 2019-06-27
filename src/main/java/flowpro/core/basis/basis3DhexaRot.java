package flowpro.core.basis;

import java.io.IOException;

/**
 *
 * @author obublik
 */
public class basis3DhexaRot extends Basis {

    public double[][] coeffs;    // koeficieny bazovych polynomu
    private int[] xExp;   //  exponent of variable xi of the k-th basis function
    private int[] yExp;  // exponent of variable eta of the k-th basis function
    private int[] zExp;  // exponent of variable zeta of the k-th basis function  
    private int order;

    public basis3DhexaRot(int order) throws IOException {
        this.order = order;
        int[] nBas = new int[]{1, 5, 10};
        nBasis = nBas[order - 1];
        //               1  2                 3
        xExp = new int[]{0, 1, 0, 0, 0, 0, 1, 2, 0, 1};
        yExp = new int[]{0, 0, 1, 0, 2, 1, 1, 0, 0, 0};
        zExp = new int[]{0, 0, 0, 1, 0, 1, 0, 0, 2, 1};
        
//        int[] nBas = new int[]{1, 4, 10, 20, 35};
//        nBasis = nBas[order - 1];
//        //               1  2        3                 4                             5
//        xExp = new int[]{0, 1, 0, 0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1};
//        yExp = new int[]{0, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1};
//        zExp = new int[]{0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2};
    }

    public void calculateCoefficients() {
        basisType = "taylor";
        coeffs = new double[nBasis][nBasis];
        for(int i = 0; i < nBasis; i++){
            coeffs[i][i] = 1;
        }
    }

    // tato funkce vraci hodnotu m-te bazove funkce v bode Xi
    public double basisFun(int m, double[] Xi) {
        double val = 0;
        for (int i = 0; i < nBasis; i++) {
            val = val + coeffs[i][m] * Math.pow(Xi[0], xExp[i]) * Math.pow(Xi[1], yExp[i]) * Math.pow(Xi[2], zExp[i]);
        }
        return val;
    }

    // tato funkce vraci hodnotu x-ove derivace m-te bazove funkce v bode x, y
    public double derBasis(int m, double[] Xi, int dim) {
        switch (dim) {
            case 0:
                double dx = 0;
                for (int i = 0; i < nBasis; i++) {
                    if (xExp[i] > 0) {
                        dx = dx + coeffs[i][m] * xExp[i] * Math.pow(Xi[0], xExp[i] - 1)
                                * Math.pow(Xi[1], yExp[i]) * Math.pow(Xi[2], zExp[i]);
                    }
                }
                return dx;
            case 1:
                double dy = 0;
                for (int i = 0; i < nBasis; i++) {
                    if (yExp[i] > 0) {
                        dy = dy + coeffs[i][m] * yExp[i] * Math.pow(Xi[0], xExp[i])
                                * Math.pow(Xi[1], yExp[i] - 1) * Math.pow(Xi[2], zExp[i]);
                    }
                }
                return dy;
            case 2:
                double dz = 0;
                for (int i = 0; i < nBasis; i++) {
                    if (zExp[i] > 0) {
                        dz = dz + coeffs[i][m] * zExp[i] * Math.pow(Xi[0], xExp[i])
                                * Math.pow(Xi[1], yExp[i]) * Math.pow(Xi[2], zExp[i] - 1);
                    }
                }
                return dz;

            default:
                return 0;
        }
    }
}
