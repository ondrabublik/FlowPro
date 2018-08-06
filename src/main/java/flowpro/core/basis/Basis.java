package flowpro.core.basis;

import java.io.Serializable;

/**
 *
 * @author ales
 */
abstract public class Basis implements Serializable{
    public int nBasis;
    public String basisType = "lagrange";

    abstract public void calculateCoefficients();
    
    abstract public double basisFun(int m, double[] Xi);
    
    abstract public double derBasis(int m, double[] Xi, int dim);
    
    public double[][] getBasisXi(double[][] coordinates) {
        int nPoints = coordinates.length;
        double[][] basis = new double[nPoints][nBasis];
        for (int p = 0; p < nPoints; p++) {
            for (int m = 0; m < nBasis; m++) {
                basis[p][m] = basisFun(m, coordinates[p]);
            }
        }
        return basis;
    }
    
    public double[][][] getDerBasisXi(double[][] coordinates, int dim) {
        int nPoints = coordinates.length;
        double[][][] derBasis = new double[nPoints][nBasis][dim];
        for (int d = 0; d < dim; d++) {
            for (int p = 0; p < nPoints; p++) {
                for (int m = 0; m < nBasis; m++) {
                    derBasis[p][m][d] = derBasis(m, coordinates[p], d);
                }
            }
        }
        return derBasis;
    }
}
