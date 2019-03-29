package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class transformation3Dhexa extends Transformation {

    private double[][] A;
    private final double[][] invV;

    public transformation3Dhexa() {
        coordsXi = new double[][]{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

        double[][] V = new double[8][8]; // Vandermondova matice
        for (int i = 0; i < 8; i++) {
            V[i][0] = coordsXi[i][0] * coordsXi[i][1] * coordsXi[i][2];
            V[i][1] = coordsXi[i][0] * coordsXi[i][1];
            V[i][2] = coordsXi[i][0] * coordsXi[i][2];
            V[i][3] = coordsXi[i][1] * coordsXi[i][2];
            V[i][4] = coordsXi[i][0];
            V[i][5] = coordsXi[i][1];
            V[i][6] = coordsXi[i][2];
            V[i][7] = 1;
        }
        invV = Mat.invert(V);
    }

    @Override
    public void computeTransform(double[][] vertices) {
        coordsX = vertices;
        A = new double[3][8];
        double[] aux = new double[vertices.length];
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = vertices[i][0];
        }
        A[0] = Mat.times(invV, aux);

        for (int i = 0; i < vertices.length; i++) {
            aux[i] = vertices[i][1];
        }
        A[1] = Mat.times(invV, aux);

        for (int i = 0; i < vertices.length; i++) {
            aux[i] = vertices[i][2];
        }
        A[2] = Mat.times(invV, aux);
    }

    // transform
    @Override
    public double[] getX(double[] Xi) {
        return new double[]{A[0][0] * Xi[0] * Xi[1] * Xi[2] + A[0][1] * Xi[0] * Xi[1] + A[0][2] * Xi[0] * Xi[2] + A[0][3] * Xi[1] * Xi[2] + A[0][4] * Xi[0] + A[0][5] * Xi[1] + A[0][6] * Xi[2] + A[0][7],
            A[1][0] * Xi[0] * Xi[1] * Xi[2] + A[1][1] * Xi[0] * Xi[1] + A[1][2] * Xi[0] * Xi[2] + A[1][3] * Xi[1] * Xi[2] + A[1][4] * Xi[0] + A[1][5] * Xi[1] + A[1][6] * Xi[2] + A[1][7],
            A[2][0] * Xi[0] * Xi[1] * Xi[2] + A[2][1] * Xi[0] * Xi[1] + A[2][2] * Xi[0] * Xi[2] + A[2][3] * Xi[1] * Xi[2] + A[2][4] * Xi[0] + A[2][5] * Xi[1] + A[2][6] * Xi[2] + A[2][7]};
    }

    
    /*
    public double[] getXi(double[] X) {
        double[] Xi = new double[]{0.5, 0.5, 0.5};
        double[] Xin = new double[]{0.5, 0.5, 0.5};
        double[] V = new double[]{0, 0, 0};
        double[][] J = new double[3][3];
        double h = 1e-10;
        for (int i = 0; i <= 500; i++) {
            for (int j = 0; j < 3; j++) {
                V[j] = h;
                J[j] = Mat.times(Mat.minusVec(getX(Mat.plusVec(Xi, V)), getX(Xi)), 1 / h);
                V[j] = 0;
            }
            Mat.transposeInPlace(J);
            Xin = Mat.minusVec(Xi, Mat.times(Mat.invert(J), Mat.minusVec(getX(Xi),X)));
            
            
            if (Mat.L2Norm(Mat.minusVec(Xin, Xi)) < 1e-11) {
                break;
            }
            if (i == 500 || Double.isNaN(Xin[0]) || Double.isNaN(Xin[1]) || Double.isNaN(Xin[2])) {
                System.out.println("Inverse mapping error!" + Mat.L2Norm(Mat.minusVec(Xin, Xi)));
                Xin = new double[]{0.5, 0.5, 0.5};
                break;
            }

            System.arraycopy(Xin, 0, Xi, 0, 3);
        }
        return Xin;
    }
    */

    @Override
    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        double D = 0;
        switch (dimBottom) {
            case 0:
                D = A[dimTop][0] * Xi[1] * Xi[2] + A[dimTop][1] * Xi[1] + A[dimTop][2] * Xi[2] + A[dimTop][4];
                break;
            case 1:
                D = A[dimTop][0] * Xi[0] * Xi[2] + A[dimTop][1] * Xi[0] + A[dimTop][3] * Xi[2] + A[dimTop][5];
                break;
            case 2:
                D = A[dimTop][0] * Xi[0] * Xi[1] + A[dimTop][2] * Xi[0] + A[dimTop][3] * Xi[1] + A[dimTop][6];
                break;
        }
        return D;
    }

    @Override
    public double jacobian(double[] Xi) {
        return Math.abs(Mat.det(jacobiMatrix(Xi)));
    }

    @Override
    public double[][] getInterpolant(Quadrature quad) {
        double[][] interpolant = new double[quad.nPoints][8];
        for (int p = 0; p < quad.nPoints; p++) {
            double xi = quad.coords[p][0];
            double eta = quad.coords[p][1];
            double zeta = quad.coords[p][2];
            interpolant[p][0] = (1 - xi) * (1 - eta) * (1 - zeta);
            interpolant[p][1] = xi * (1 - eta) * (1 - zeta);
            interpolant[p][2] = xi * eta * (1 - zeta);
            interpolant[p][3] = (1 - xi) * eta * (1 - zeta);
            interpolant[p][4] = (1 - xi) * (1 - eta) * zeta;
            interpolant[p][5] = xi * (1 - eta) * zeta;
            interpolant[p][6] = xi * eta * zeta;
            interpolant[p][7] = (1 - xi) * eta * zeta;
        }
        return interpolant;
    }
    
    @Override
    public double[] getXs(){
        return getX(new double[]{0.5, 0.5, 0.5});
    }
}
