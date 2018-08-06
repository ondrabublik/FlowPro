package flowpro.core.transformation;

import flowpro.core.quadrature.Quadrature;

public class faceTransformation1Dpoint extends FaceTransformation {

    // public static final int nk = 3;
    private static final double[] coordsXi = {1};
    private double[][] coordsXiRef;
    
    public faceTransformation1Dpoint(Transformation transform, int[] faceIndexes) {
        this.transform = transform;
        coordsXiRef = new double[1][1];
        coordsXiRef[0][0] = transform.coordsXi[faceIndexes[0]][0];
    }
    
    @Override
    public double[] getXiRef(double[] Xi) {
        return new double[]{coordsXiRef[0][0]};
    }

    @Override
    public double[] getX(double[] Xi) {
        return transform.getX(getXiRef(Xi));
    }

    @Override
    public double jacobian(double[] Xi) {
        return 1.0;
    }

    @Override
    public double[][] getInterpolant(Quadrature quad) {
        return new double[][]{{1}};
    }

    @Override
    public double[] getNormal(double[] Xi) {
        return new double[]{1};
    }
}

