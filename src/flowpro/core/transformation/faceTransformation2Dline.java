package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

public class faceTransformation2Dline extends FaceTransformation {

    // public static final int nk = 3;
    private static final double[] coordsXi = {-1, 1};
    private double[][] coordsXiRef;

    public faceTransformation2Dline(Transformation transform, int[] faceIndexes) {
        this.transform = transform;
        coordsXiRef = new double[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                coordsXiRef[i][j] = transform.coordsXi[faceIndexes[i]][j];
            }
        }
    }

    @Override
    public double[] getXiRef(double[] Xi) {
        double[] XiRef = new double[2];
        XiRef[0] = (1 + Xi[0]) / 2 * coordsXiRef[1][0] + (1 - Xi[0]) / 2 * coordsXiRef[0][0];
        XiRef[1] = (1 + Xi[0]) / 2 * coordsXiRef[1][1] + (1 - Xi[0]) / 2 * coordsXiRef[0][1];
        return XiRef;
    }

    @Override
    public double[] getX(double[] Xi) {
        return transform.getX(getXiRef(Xi));
    }

    @Override
    public double jacobian(double[] Xi) {
        double[] s = new double[]{coordsXiRef[1][0] - coordsXiRef[0][0], coordsXiRef[1][1] - coordsXiRef[0][1]};
        double[] Js = transform.jacobianMatrixDirection(getXiRef(Xi), s);
        return Mat.L2Norm(Js) / 2;
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

    @Override
    public double[] getNormal(double[] Xi) {
        double[] s = new double[]{coordsXiRef[1][0] - coordsXiRef[0][0], coordsXiRef[1][1] - coordsXiRef[0][1]};
        double[] Js = transform.jacobianMatrixDirection(getXiRef(Xi), s);
        double[] n = Mat.normVector(new double[]{Js[1], -Js[0]});
        return n;
    }
}
