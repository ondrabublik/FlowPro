package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class transformation3Dtetra extends Transformation {

    // public static final int nk = 3;
    private double[] A, B, C;
    private double[] iA, iB, iC;
    private double[][] invV, invVi;

    public transformation3Dtetra(double[][] vertices) {
        coordsXi = new double[][]{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        coordsX = vertices;
        double[][] V;
        V = new double[4][4];
        for (int i = 0; i < 4; i++) {
            V[i][0] = coordsXi[i][0];
            V[i][1] = coordsXi[i][1];
            V[i][2] = coordsXi[i][2];
            V[i][3] = 1;
        }
        invV = Mat.invert(V);

        for (int i = 0; i < 4; i++) {
            V[i][0] = coordsX[i][0];
            V[i][1] = coordsX[i][1];
            V[i][2] = coordsX[i][2];
            V[i][3] = 1;
        }
        invVi = Mat.invert(V);
    }

    @Override
    public void computeTransform(double[][] vertices) {
        coordsX = vertices;
        double[] aux = new double[vertices.length];
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = coordsX[i][0];
        }
        A = Mat.times(invV, aux);
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = coordsX[i][1];
        }
        B = Mat.times(invV, aux);
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = coordsX[i][2];
        }
        C = Mat.times(invV, aux);

        // inverse transform
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = coordsXi[i][0];
        }
        iA = Mat.times(invVi, aux);
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = coordsXi[i][1];
        }
        iB = Mat.times(invVi, aux);
        for (int i = 0; i < vertices.length; i++) {
            aux[i] = coordsXi[i][2];
        }
        iC = Mat.times(invVi, aux);
    }

    @Override
    public double[] getX(double[] Xi) {
        return new double[]{A[0] * Xi[0] + A[1] * Xi[1] + A[2] * Xi[2] + A[3], B[0] * Xi[0] + B[1] * Xi[1] + B[2] * Xi[2] + B[3], C[0] * Xi[0] + C[1] * Xi[1] + C[2] * Xi[2] + C[3]};
    }

    @Override
    public double[] getXi(double[] X) {
        return new double[]{iA[0] * X[0] + iA[1] * X[1] + iA[2] * X[2] + iA[3], iB[0] * X[0] + iB[1] * X[1] + iB[2] * X[2] + iB[3], iC[0] * X[0] + iC[1] * X[1] + iC[2] * X[2] + iC[3]};
    }

    @Override
    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        double D = 0;
        switch (dimTop) {
            case 0:
                D = A[dimBottom];
                break;
            case 1:
                D = B[dimBottom];
                break;
            case 2:
                D = C[dimBottom];
                break;
        }
        return D;
    }

    @Override
    public double jacobian(double[] Xi) {
        return Math.abs(A[0] * B[1] * C[2] + A[1] * B[2] * C[0] + A[2] * B[0] * C[1] - A[2] * B[1] * C[0] - A[1] * B[0] * C[2] - A[0] * B[2] * C[1]);
    }

    @Override
    public double[][] getInterpolant(Quadrature quad) {
        double[][] interpolant = new double[quad.nPoints][4];
        for (int p = 0; p < quad.nPoints; p++) {
            double xi = quad.coords[p][0];
            double eta = quad.coords[p][1];
            double zeta = quad.coords[p][2];
            interpolant[p][0] = 1 - xi - eta - zeta;
            interpolant[p][1] = xi;
            interpolant[p][2] = eta;
            interpolant[p][3] = zeta;
        }
        return interpolant;
    }
}
