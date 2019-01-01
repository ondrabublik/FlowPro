package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;
import flowpro.core.curvedBoundary.Curved2DLine;
import flowpro.core.curvedBoundary.FaceCurvature;

public class transformation2DTriangleCurved extends Transformation {

    Curved2DLine fCurv;
    
    // public static final int nk = 3;
    private double[] A, B;
    private double[][] invV;

    public transformation2DTriangleCurved(FaceCurvature fCurv) {
        this.fCurv = (Curved2DLine)fCurv;
        coordsXi = new double[][]{{0, 1}, {0, 0}, {1, 0}, {0.5, 0.5}};

        // transform
        double[][] V;
        V = new double[4][4];
        for (int i = 0; i < 4; i++) {
            V[i][0] = coordsXi[i][0] * coordsXi[i][1];
            V[i][1] = coordsXi[i][0];
            V[i][2] = coordsXi[i][1];
            V[i][3] = 1;
        }
        invV = Mat.invert(V);
    }

    @Override
    public void computeTransform(double[][] vertices) {
        coordsX = new double[vertices.length + 1][2];
        //double[] Xr = verticeRadius(vertices[0],vertices[2],1.3);
        double[] Xr = fCurv.getPoint(0.5);
        for (int i = 0; i < 2; i++) {
            coordsX[0][i] = vertices[0][i];
            coordsX[1][i] = vertices[1][i];
            coordsX[2][i] = vertices[2][i];
            coordsX[3][i] = Xr[i];
        }
        double[] aux = new double[coordsX.length];
        for (int i = 0; i < coordsX.length; i++) {
            aux[i] = coordsX[i][0];
        }
        A = Mat.times(invV, aux);
        for (int i = 0; i < coordsX.length; i++) {
            aux[i] = coordsX[i][1];
        }
        B = Mat.times(invV, aux);
    }

    @Override
    public double[] getX(double[] Xi) {
        return new double[]{A[0] * Xi[0]* Xi[1] + A[1] * Xi[0] + A[2]* Xi[1] + A[3], B[0] * Xi[0]* Xi[1] + B[1] * Xi[0] + B[2]* Xi[1] + B[3]};
    }

    @Override
    public double[] getXi(double[] X) {
        double[] Xi = new double[]{0.25, 0.25};
        double[] Xin = new double[]{0.25, 0.25};
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
            Xin = Mat.minusVec(Xi, Mat.times(Mat.invert(J), Mat.minusVec(getX(Xi), X)));

            if (Mat.L2Norm(Mat.minusVec(Xin, Xi)) < 1e-13) {
                break;
            }

            if (i == 50 || Double.isNaN(Xin[0]) || Double.isNaN(Xin[1])) {
                System.out.println("Inverse mapping error!");
                Xin = new double[]{0.25, 0.25};
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
    
    double[] verticeRadius(double[] x1, double[] x2, double r){
        double[] xs = Mat.times(Mat.plusVec(x1,x2),0.5);
        double L = Mat.L2Norm(Mat.minusVec(x1,x2));
        double h = r-Math.pow(r*r-L*L/4,0.5);
        double[] n = new double[]{x1[1]-x2[1],x2[0]-x1[0]};
        n = Mat.times(n, 1/Mat.L2Norm(n));
        return Mat.plusVec(xs, Mat.times(n,-h));
    }
    
    @Override
    public double[] getXs(){
        return getX(new double[]{1.0/3,1.0/3});
    }
}