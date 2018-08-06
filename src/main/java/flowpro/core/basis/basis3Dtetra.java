package flowpro.core.basis;

import flowpro.api.Mat;
import java.io.IOException;

/**
 *
 * @author ales
 */
public class basis3Dtetra extends Basis {

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
    public basis3Dtetra(int order) throws IOException {
        this.order = order;
        nBasis = (3 * (order - 1) * order + 2) / 2;
        //               1  2        3                 4
        xExp = new int[]{0, 1, 0, 0, 1, 1, 0, 2, 0, 0, 2, 2, 1, 0, 1, 0, 3, 0, 0};
        yExp = new int[]{0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 2, 0, 1, 0, 3, 0};
        zExp = new int[]{0, 0, 0, 1, 0, 1, 1, 0, 0, 2, 0, 1, 0, 1, 2, 2, 0, 0, 3};
    }

    public void calculateCoefficients() {
        if (order == 1) {
            double[][] coefficients = new double[1][1];
            coefficients[0][0] = 1;
            coeffs = coefficients;
        } else {
            double[][] b = barycentr(order, nBasis);
            double[][] V = new double[nBasis][nBasis];
            for (int i = 0; i < nBasis; i++) {
                double[] X = new double[]{b[i][1], b[i][2], b[i][3]};
                for (int j = 0; j < nBasis; j++) {
                    V[i][j] = Math.pow(X[0], xExp[j]) * Math.pow(X[1], yExp[j]) * Math.pow(X[2], zExp[j]);
                }
            }
            coeffs = Mat.invert(V);
        }
    }

    double[][] barycentr(int rad, int ne) {
        double[][] b;
        if (rad == 1) {
            b = new double[1][4];
            b[0][0] = 1.0 / 4;
            b[0][1] = 1.0 / 4;
            b[0][2] = 1.0 / 4;
            b[0][3] = 1.0 / 4;
        } else {
            b = new double[][]{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        }
        return b;
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
