package flowpro.core.basis;

import flowpro.api.Mat;
import java.io.IOException;

/**
 *
 * @author ales
 */
public class basis1Dline extends Basis {

    public double[][] coeffs;    // koeficieny bazovych polynomu
    private int[] xExp;   //  exponent of variable xi of the k-th basis function
    private int order;

    /**
     *
     * @param order
     * @throws IOException
     */
    public basis1Dline(int order) throws IOException {
        this.order = order;
        nBasis = order;

        xExp = new int[order];
        for (int i = 0; i < order; i++) {
            xExp[i] = i;
        }
    }

    public void calculateCoefficients() {
        if (order == 1) {
            double[][] coefficients = new double[1][1];
            coefficients[0][0] = 1;
            coeffs = coefficients;
        } else {
            double[][] V = new double[nBasis][nBasis];
            for (int i = 0; i < nBasis; i++) {
                double[] X = new double[]{-1.0 + (2.0*i)/(nBasis-1)};
                for (int j = 0; j < nBasis; j++) {
                    V[i][j] = Math.pow(X[0], xExp[j]);
                }
            }
            coeffs = Mat.invert(V);
        }
    }

    // tato funkce vraci hodnotu m-te bazove funkce v bode X
    public double basisFun(int m, double[] X) {
        double val = 0;
        for (int i = 0; i < nBasis; i++) {
            val = val + coeffs[i][m] * Math.pow(X[0], xExp[i]);
        }
        return val;
    }

    // tato funkce vraci hodnotu x-ove derivace m-te bazove funkce v bode x, y
    public double derBasis(int m, double[] X, int dim) {
        double dx = 0;
        for (int i = 0; i < nBasis; i++) {
            if (xExp[i] > 0) {
                dx = dx + coeffs[i][m] * xExp[i] * Math.pow(X[0], xExp[i] - 1);
            }
        }
        return dx;
    }
}
