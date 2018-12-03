package flowpro.core.basis;

import flowpro.api.Mat;
import java.io.IOException;

/**
 *
 * @author ales
 */
public class basis2DsquareTaylor extends Basis {

    public double[][] coeffs;    // koeficieny bazovych polynomu
    private final int[] xExp;   //  exponent of variable xi of the k-th basis function
    private final int[] yExp;  // exponent of variable eta of the k-th basis function   
    private final int order;

    public basis2DsquareTaylor(int order) throws IOException {
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

     public void calculateCoefficients() {
        basisType = "taylor";
        if (order == 1) {
            double[][] coefficients = new double[1][1];
            coefficients[0][0] = 1;
            coeffs = coefficients;
        } else {
            coeffs = new double[nBasis][nBasis];
            for(int i = 0; i < nBasis; i++){
                coeffs[i][i] = 1;
            }
        }
    }
    
    // tato funkce vraci hodnotu m-te bazove funkce v bode Xi
    public double basisFun(int m, double[] Xi) {
        double val = 0;
        for (int i = 0; i < nBasis; i++) {
            val = val + coeffs[i][m] * Math.pow(Xi[0], xExp[i]) * Math.pow(Xi[1], yExp[i]);
        }
        return val;
    }

    // tato funkce vraci hodnotu x-ove derivace m-te bazove funkce v bode x, y
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

