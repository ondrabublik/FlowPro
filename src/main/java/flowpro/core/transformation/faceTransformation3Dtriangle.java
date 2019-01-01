package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;

/**
 *
 * @author obublik
 */
public class faceTransformation3Dtriangle extends FaceTransformation{

    // public static final int nk = 3;
    private double[][] coordsXi = new double[][]{{0,0},{1,0},{0,1}};
    private double[][] coordsXiRef;
    private double[] A, B, C;
    private double[][] invV;

    public faceTransformation3Dtriangle(Transformation transform, int[] faceIndexes) {
        this.transform = transform;
        coordsXiRef = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                coordsXiRef[i][j] = transform.coordsXi[faceIndexes[i]][j];
            }
        }
        
        double[][] V;
        V = new double[3][3];
        for (int i = 0; i < 3; i++) {
            V[i][0] = coordsXi[i][0];
            V[i][1] = coordsXi[i][1];
            V[i][2] = 1;
        }
        invV = Mat.invert(V);
    
        double[] aux = new double[coordsXiRef.length]; 
        for (int i = 0; i < coordsXiRef.length; i++) {
            aux[i] = coordsXiRef[i][0];
        }
        A = Mat.times(invV, aux);
        for (int i = 0; i < coordsXiRef.length; i++) {
            aux[i] = coordsXiRef[i][1];
        }
        B = Mat.times(invV, aux);
        for (int i = 0; i < coordsXiRef.length; i++) {
            aux[i] = coordsXiRef[i][2];
        }
        C = Mat.times(invV, aux);
    }

    @Override
    public double[] getXiRef(double[] Xi) {
        return new double[]{A[0] * Xi[0] + A[1] * Xi[1] + A[2], B[0] * Xi[0] + B[1] * Xi[1] + B[2], C[0] * Xi[0] + C[1] * Xi[1] + C[2]};
    }

    @Override
    public double[] getX(double[] Xi) {
        return transform.getX(getXiRef(Xi));
    }

    @Override
    public double jacobian(double[] Xi) {
        double[] s1 = new double[]{A[0], B[0], C[0]};
        double[] s2 = new double[]{A[1], B[1], C[1]};
        double[] Js1 = transform.jacobianMatrixDirection(getXiRef(Xi), s1);
        double[] Js2 = transform.jacobianMatrixDirection(getXiRef(Xi), s2);
        return Mat.L2Norm(Mat.cross(Js1,Js2))/2;
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
    public double[] getNormal(double[] Xi) {
        double[] s1 = new double[]{A[0], B[0], C[0]};
        double[] s2 = new double[]{A[1], B[1], C[1]};
        double[] Js1 = transform.jacobianMatrixDirection(getXiRef(Xi), s1);
        double[] Js2 = transform.jacobianMatrixDirection(getXiRef(Xi), s2);
        return Mat.normVector(Mat.cross(Js2,Js1));
    }
    
    @Override
    public double[] getXs(){
        return getX(new double[]{1.0/3, 1.0/3});
    }
}

