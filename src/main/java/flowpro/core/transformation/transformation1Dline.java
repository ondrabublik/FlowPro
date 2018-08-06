package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class transformation1Dline extends Transformation {

    // public static final int nk = 3;
    private double[] A;
    private final double[][] invV;

    public transformation1Dline() {
        coordsXi = new double[2][1];
        coordsXi[0][0] = -1;
        coordsXi[1][0] = 1;
        
        double[][] V;
        V = new double[2][2];
        for (int i = 0; i < 2; i++) {
            V[i][0] = coordsXi[i][0];
            V[i][1] = 1;
        }
        invV = Mat.invert(V);
    }

    @Override
    public void computeTransform(double[][] vertices){
        coordsX = vertices;
        double[] aux = new double[vertices.length];
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = vertices[i][0];
        }
        A = Mat.times(invV, aux);
    }
    
    @Override
    public double[] getX(double[] Xi) {
        return new double[]{A[0] * Xi[0] + A[1]};
    }

    @Override
    public double[] getXi(double[] X) {
        return new double[]{(X[0]-A[1])/A[0]};
    }
    
    @Override
    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        return A[0];
    }
    
    @Override
    public double jacobian(double Xi[]) {
        return Math.abs(A[0]);
    }
    
    @Override
    public double[][] getInterpolant(Quadrature quad) {
        double[][] interpolant = new double[quad.nPoints][2];

        for (int p = 0; p < quad.nPoints; p++) {
            double xi = quad.coords[p][0];
            interpolant[p][0] = (1 - xi) / 2;
            interpolant[p][1] = (1 + xi) / 2;
        }
        return interpolant;
    }
}

