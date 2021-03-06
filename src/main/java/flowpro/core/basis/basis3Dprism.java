
package flowpro.core.basis;

import flowpro.api.Mat;
import java.io.IOException;


public class basis3Dprism extends Basis {

    public double[][] coeffs;    // koeficieny bazovych polynomu
    private int[] xExp;   //  exponent of variable xi of the k-th basis function
    private int[] yExp;  // exponent of variable eta of the k-th basis function
    private int[] zExp;  // exponent of variable eta of the k-th basis function
    private int order;

    /**
     *
     * @param order
     * @throws IOException
     */
    public basis3Dprism(int order) throws IOException {
        this.order = order;
        int[] nBas = new int[]{1, 4};
        nBasis = nBas[order - 1];

        xExp = new int[]{0, 1, 0, 0, 1, 1, 0, 2, 0, 0};
        yExp = new int[]{0, 0, 1, 0, 1, 0, 1, 0, 2, 0};
        zExp = new int[]{0, 0, 0, 1, 0, 1, 1, 0, 0, 2};
    }

    public void calculateCoefficients() {
        if (order == 1) {
            double[][] coefficients = new double[1][1];
            coefficients[0][0] = 1;
            coeffs = coefficients;
        } else {
            double[] xi =   new double[]{1./3, 0.5,   0, 1./3};
            double[] eta =  new double[]{1./3,   0, 0.5, 1./3};
            double[] zeta = new double[]{ 0.5,   0,   0,    1};
            double[][] V = new double[nBasis][nBasis];
            for (int i = 0; i < nBasis; i++) {
                for (int j = 0; j < nBasis; j++) {
                    V[i][j] = Math.pow(xi[i], xExp[j]) * Math.pow(eta[i], yExp[j]) * Math.pow(zeta[i], zExp[j]);
                }
            }
            coeffs = Mat.invert(V);
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