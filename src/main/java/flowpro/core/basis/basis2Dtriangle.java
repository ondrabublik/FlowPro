package flowpro.core.basis;

import flowpro.api.Mat;
import java.io.IOException;

/**
 *
 * @author ales
 */
public class basis2Dtriangle extends Basis {

    public double[][] coeffs;    // koeficieny bazovych polynomu
    public int[] xExp;   //  exponent of variable xi of the k-th basis function
    public int[] yExp;  // exponent of variable eta of the k-th basis function
    public int order;

    /**
     *
     * @param order
     * @throws IOException
     */
    public basis2Dtriangle(int order) throws IOException {
        this.order = order;
        nBasis = order * (1 + order) / 2;

        xExp = new int[order * order];
        yExp = new int[order * order];
        int s = 0;
        for (int i = 0; i < order; i++) {
            for (int j = 0; j <= i; j++) {
                xExp[s] = i - j;
                yExp[s] = j;
                ++s;
            }
        }
        for (int i = order - 1; i > 0; i--) {
            for (int j = 0; j < i; j++) {
                xExp[s] = order - 1 - j;
                yExp[s] = order - i + j;
                ++s;
            }
        }
    }

    @Override
    public void calculateCoefficients() {
        if (order == 1) {
            double[][] coefficients = new double[1][1];
            coefficients[0][0] = 1;
            coeffs = coefficients;
        } else {
            double[][] b = barycentr(order, nBasis);
            double[][] V = new double[nBasis][nBasis];
            for (int i = 0; i < nBasis; i++) {
                double[] X = new double[]{b[i][1], b[i][2]};
                for (int j = 0; j < nBasis; j++) {
                    V[i][j] = Math.pow(X[0], xExp[j]) * Math.pow(X[1], yExp[j]);
                }
            }
            coeffs = Mat.invert(V);
        }
    }

    double[][] barycentr(int rad, int ne) {
        double[][] b;
        if (rad == 1) {
            b = new double[1][3];
            b[0][0] = 1.0 / 3;
            b[0][1] = 1.0 / 3;
            b[0][2] = 1.0 / 3;
        } else {
            double[] x = Mat.linspace(0, 1, rad);
            b = new double[ne][3];
            int s = 0;
            for (int i = 0; i < rad; i++) {
                double[] y = Mat.linspace(0, x[i], i + 1);
                for (int j = 0; j < i + 1; j++) {
                    b[s][0] = 1 - x[i];
                    b[s][1] = y[j];
                    b[s][2] = x[i] - y[j];
                    s = s + 1;
                }
            }
        }
        return b;
    }

    // tato funkce vraci hodnotu m-te bazove funkce v bode Xi
    @Override
    public double basisFun(int m, double[] Xi) {
        double val = 0;
        for (int i = 0; i < nBasis; i++) {
            val = val + coeffs[i][m] * Math.pow(Xi[0], xExp[i]) * Math.pow(Xi[1], yExp[i]);
        }
        return val;
    }

    // tato funkce vraci hodnotu x-ove derivace m-te bazove funkce v bode x, y
    @Override
    public double derBasis(int m, double[] Xi, int dim) {
        double der = 0;
        if (dim == 0) {
            for (int i = 0; i < nBasis; i++) {
                if (xExp[i] > 0) {
                    der = der + coeffs[i][m] * xExp[i] * Math.pow(Xi[0], xExp[i] - 1)
                            * Math.pow(Xi[1], yExp[i]);
                }
            }
        } else {
            for (int i = 0; i < nBasis; i++) {
                if (yExp[i] > 0) {
                    der = der + coeffs[i][m] * yExp[i] * Math.pow(Xi[0], xExp[i])
                            * Math.pow(Xi[1], yExp[i] - 1);
                }
            }
        }
        return der;
    }
}
