package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class transformation3DhexaCylinderCoords extends Transformation {

    private double[][] A;
    private final double[][] invV;

    public transformation3DhexaCylinderCoords() {
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
        double[] Xr = new double[]{A[0][0] * Xi[0] * Xi[1] * Xi[2] + A[0][1] * Xi[0] * Xi[1] + A[0][2] * Xi[0] * Xi[2] + A[0][3] * Xi[1] * Xi[2] + A[0][4] * Xi[0] + A[0][5] * Xi[1] + A[0][6] * Xi[2] + A[0][7],
            A[1][0] * Xi[0] * Xi[1] * Xi[2] + A[1][1] * Xi[0] * Xi[1] + A[1][2] * Xi[0] * Xi[2] + A[1][3] * Xi[1] * Xi[2] + A[1][4] * Xi[0] + A[1][5] * Xi[1] + A[1][6] * Xi[2] + A[1][7],
            A[2][0] * Xi[0] * Xi[1] * Xi[2] + A[2][1] * Xi[0] * Xi[1] + A[2][2] * Xi[0] * Xi[2] + A[2][3] * Xi[1] * Xi[2] + A[2][4] * Xi[0] + A[2][5] * Xi[1] + A[2][6] * Xi[2] + A[2][7]};

        return cylinderTransform(Xr);
    }

    @Override
    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        double[] Y = getX(Xi);
        double r = Y[2];
        double fi = Y[1];
        double dx = A[dimTop][0] * Xi[1] * Xi[2] + A[dimTop][1] * Xi[1] + A[dimTop][2] * Xi[2] + A[dimTop][4];
        double dr = A[dimTop][0] * Xi[0] * Xi[2] + A[dimTop][1] * Xi[0] + A[dimTop][3] * Xi[2] + A[dimTop][5];
        double dfi = A[dimTop][0] * Xi[0] * Xi[1] + A[dimTop][2] * Xi[0] + A[dimTop][3] * Xi[1] + A[dimTop][6];

        double D = 0;
        switch (dimBottom) {
            case 0:
                D = dx;
                break;
            case 1:
                D = dr * Math.cos(fi) - r * Math.sin(fi) * dfi;
                break;
            case 2:
                D = -dr * Math.sin(fi) - r * Math.cos(fi) * dfi;
                break;
        }
        return D;
    }

    public double[] cylinderTransform(double[] X) {
        return new double[]{X[0], X[2] * Math.cos(X[1]), -X[2] * Math.sin(X[1])};
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
