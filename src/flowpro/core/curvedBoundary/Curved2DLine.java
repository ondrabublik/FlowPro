package flowpro.core.curvedBoundary;

import flowpro.api.Mat;

/**
 *
 * @author obublik
 */
public class Curved2DLine extends FaceCurvature {

    double[] A;
    double[] B;

    Curved2DLine(double[] t0, double[] t1, double[] X0, double[] X1) {
        double[][] V = new double[][]{{0, 0, 0, 1}, {1, 1, 1, 1}, {0, 0, 1, 0}, {3, 2, 1, 0}};
        double[][] iV = Mat.invert(V);

        double L = Mat.L2Norm(Mat.minusVec(X0, X1));
        double[] aux = new double[]{X0[0], X1[0], t0[0] * L, t1[0] * L};
        A = Mat.times(iV, aux);

        aux = new double[]{X0[1], X1[1], t0[1] * L, t1[1] * L};
        B = Mat.times(iV, aux);
    }

    Curved2DLine(String side, double[] t, double[] X0, double[] X1) {
        double[][] V;
        if ("left".equals(side)) {
            V = new double[][]{{0, 0, 1}, {1, 1, 1}, {0, 1, 0}};
        } else {
            V = new double[][]{{0, 0, 1}, {1, 1, 1}, {2, 1, 0}};
        }
        double[][] iV = Mat.invert(V);

        double L = Mat.L2Norm(Mat.minusVec(X0, X1));
        double[] aux = new double[]{X0[0], X1[0], t[0] * L};
        A = Mat.times(iV, aux);
        A = new double[]{0, A[0], A[1], A[2]};

        aux = new double[]{X0[1], X1[1], t[1] * L};
        B = Mat.times(iV, aux);
        B = new double[]{0, B[0], B[1], B[2]};
    }

    public double[] getPoint(double t) {
        return new double[]{A[0] * t * t * t + A[1] * t * t + A[2] * t + A[3], B[0] * t * t * t + B[1] * t * t + B[2] * t + B[3]};
    }
}
