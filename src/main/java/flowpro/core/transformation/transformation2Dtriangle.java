package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class transformation2Dtriangle extends Transformation {

    // public static final int nk = 3;
    private double[] A, B;
    private double[] iA, iB;
    private double[][] invV, invVi;

    public transformation2Dtriangle(double[][] vertices) {
        coordsXi = new double[][]{{0, 0}, {1, 0}, {0, 1}};
        coordsX = vertices;
        
        // transform
        double[][] V;
        V = new double[3][3];
        for (int i = 0; i < 3; i++) {
            V[i][0] = coordsXi[i][0];
            V[i][1] = coordsXi[i][1];
            V[i][2] = 1;
        }
        invV = Mat.invert(V);
        
        // inverse transform
        double[][] Vi = new double[3][3];
        for (int i = 0; i < 3; i++) {
            Vi[i][0] = coordsX[i][0];
            Vi[i][1] = coordsX[i][1];
            Vi[i][2] = 1;
        }
        invVi = Mat.invert(Vi);
    }

    @Override
    public void computeTransform(double[][] vertices){
        coordsX = vertices;
        double[] aux = new double[coordsX.length];
        for (int i = 0; i < coordsX.length; i++) {
            aux[i] = coordsX[i][0];
        }
        A = Mat.times(invV, aux);
        for (int i = 0; i < coordsX.length; i++) {
            aux[i] = coordsX[i][1];
        }
        B = Mat.times(invV, aux);
        
        aux = new double[coordsXi.length];
        for (int i = 0; i < coordsXi.length; i++) {
            aux[i] = coordsXi[i][0];
        }
        iA = Mat.times(invVi, aux);
        for (int i = 0; i < coordsXi.length; i++) {
            aux[i] = coordsXi[i][1];
        }
        iB = Mat.times(invVi, aux);
    }
    
    @Override
    public double[] getX(double[] Xi) {
        return new double[]{A[0] * Xi[0] + A[1] * Xi[1] + A[2], B[0] * Xi[0] + B[1] * Xi[1] + B[2]};
    }

    @Override
    public double[] getXi(double[] X) {
        return new double[]{iA[0] * X[0] + iA[1] * X[1] + iA[2], iB[0] * X[0] + iB[1] * X[1] + iB[2]};
    }

    @Override
    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        double D;
        if (dimTop == 0) {
            if (dimBottom == 0) {
                D = A[0];
            } else {
                D = A[1];
            }
        } else {
            if (dimBottom == 0) {
                D = B[0];
            } else {
                D = B[1];
            }
        }
       
        return D;
    }

    @Override
    public double jacobian(double[] Xi) {
        double dxdxi = getDX(Xi, 0, 0);
        double dxdeta = getDX(Xi, 0, 1);
        double dydxi = getDX(Xi, 1, 0);
        double dydeta = getDX(Xi, 1, 1);
        return Math.abs(dxdxi * dydeta - dxdeta * dydxi) / 2;
    }

    @Override
    public double[][] getInterpolant(Quadrature quad) {
        double[][] interpolant = new double[quad.nPoints][3];
        for (int p = 0; p < quad.nPoints; p++) {
            double xi = quad.coords[p][0];
            double eta = quad.coords[p][1];
            interpolant[p][0] = 1 - xi - eta;
            interpolant[p][1] = xi;
            interpolant[p][2] = eta;

        }
        return interpolant;
    }
    
    @Override
    public double[] getXs(){
        return getX(new double[]{1.0/3,1.0/3});
    }
}
