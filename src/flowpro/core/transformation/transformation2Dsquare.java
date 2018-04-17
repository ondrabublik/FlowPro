package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class transformation2Dsquare extends Transformation {

    private double[] A, B;
    private final double[][] invV;

    public transformation2Dsquare() {
        coordsXi = new double[][]{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        
        double[][] V = new double[4][4]; // Vandermondova matice
        for (int i = 0; i < 4; i++) {
            V[i][0] = coordsXi[i][0] * coordsXi[i][1];
            V[i][1] = coordsXi[i][0];
            V[i][2] = coordsXi[i][1];
            V[i][3] = 1;
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
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = vertices[i][1];
        }
        B = Mat.times(invV, aux);
    }
    
    // transform
    @Override
    public double[] getX(double[] Xi) {
        return new double[]{A[0] * Xi[0] * Xi[1] + A[1] * Xi[0] + A[2] * Xi[1] + A[3], B[0] * Xi[0] * Xi[1] + B[1] * Xi[0] + B[2] * Xi[1] + B[3]};
    }

    /*public double[] getXi(double[] X) {
        double xi = 0.5;
        double eta = 0.5;
        for (int i = 0; i < 50; i++) {
            double xin = (X[0] - A[3] - A[2] * eta) / (A[1] + A[0] * eta);
            double etan = (X[1] - B[3] - B[1] * xi) / (B[2] + B[0] * xi);
            double rezid = Math.abs(xin - xi) + Math.abs(etan - eta);
            xi = xin;
            eta = etan;
            if (rezid < 1e-13) {
                break;
            }
            if (i == 50 || Double.isNaN(xin) || Double.isNaN(etan)) {
                System.out.println("Inverse mapping error!");
            }
        }
        return new double[]{xi, eta};
    }*/
    
    public double[] getXi(double[] X) {
        double[] Xi = new double[]{0.5, 0.5};
        double[] Xin = new double[]{0.5, 0.5};
        double[] V = new double[]{0, 0};
        double[][] J = new double[2][2];
        double h = 1e-8;
        for (int i = 0; i < 50; i++) {
            for (int j = 0; j < 2; j++) {
                V[j] = h;
                J[j] = Mat.times(Mat.minusVec(getX(Mat.plusVec(Xi, V)), getX(Xi)), 1 / h);
                V[j] = 0;
            }
            Mat.transposeInPlace(J);
            Xin = Mat.minusVec(Xi, Mat.times(Mat.invert(J), Mat.minusVec(getX(Xi),X)));
            
            
            if (Mat.L2Norm(Mat.minusVec(Xin, Xi)) < 1e-13) {
                break;
            }
            if (i == 50 || Double.isNaN(Xin[0]) || Double.isNaN(Xin[1])) {
                System.out.println("Inverse mapping error!");
                Xin = new double[]{0.5, 0.5};
                break;
            }

            System.arraycopy(Xin, 0, Xi, 0, 2);
        }
        return Xin;
    }

    @Override
    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        double D = 0;
        if (dimTop == 0) {
            if (dimBottom == 0) {
                D = A[0] * Xi[1] + A[1];
            } else {
                D = A[0] * Xi[0] + A[2];
            }
        } else {
            if (dimBottom == 0) {
                D = B[0] * Xi[1] + B[1];
            } else {
                D = B[0] * Xi[0] + B[2];
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
        return Math.abs(dxdxi * dydeta - dxdeta * dydxi);
    }

    @Override
    public double[][] getInterpolant(Quadrature quad) {
        double[][] interpolant = new double[quad.nPoints][4];
        for (int p = 0; p < quad.nPoints; p++) {
            double xi = quad.coords[p][0];
            double eta = quad.coords[p][1];
            interpolant[p][0] = (1 - xi) * (1 - eta);
            interpolant[p][1] = xi * (1 - eta);
            interpolant[p][2] = xi * eta;
            interpolant[p][3] = (1 - xi) * eta;
        }
        return interpolant;
    }
}
