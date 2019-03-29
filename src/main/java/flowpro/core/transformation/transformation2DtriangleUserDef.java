package flowpro.core.transformation;

import flowpro.api.Complex;
import static flowpro.api.Complex.multiply;
import flowpro.api.DomainTransformationObject;
import flowpro.api.Mat;
import flowpro.core.Parameters;
import flowpro.core.quadrature.Quadrature;

public class transformation2DtriangleUserDef extends Transformation {

    // public static final int nk = 3;
    private double[] A, B;
    private final double[][] invV;
    DomainTransformationObject domTrans;

    public transformation2DtriangleUserDef(double[][] vertices, Parameters par) {
        domTrans = par.domainTransformationObject;
        if(domTrans == null){
            throw new UnsupportedOperationException("Domain transformation object not defined!");
        }
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
    }

    @Override
    public void computeTransform(double[][] vertices) {
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
    }

    @Override
    public double[] getX(double[] Xi) {
        double[] Xr = new double[]{A[0] * Xi[0] + A[1] * Xi[1] + A[2], B[0] * Xi[0] + B[1] * Xi[1] + B[2]};
        return domTrans.transform(Xr);
    }
    
    public Complex[] getX(Complex[] Xi) {
        Complex[] Xr = new Complex[Xi.length];
        Complex A0 = new Complex(A[0], 0);
        Complex A1 = new Complex(A[1], 0);
        Complex A2 = new Complex(A[2], 0);
        Xr[0] = new Complex();
        Xr[0].add(multiply(Xi[0],A0));
        Xr[0].add(multiply(Xi[1],A1));
        Xr[0].add(A2);
        
        Complex B0 = new Complex(B[0],0);
        Complex B1 = new Complex(B[1],0);
        Complex B2 = new Complex(B[2],0);
        Xr[1] = new Complex();
        Xr[1].add(multiply(Xi[0], B0));
        Xr[1].add(multiply(Xi[1], B1));
        Xr[1].add(B2);

        return domTrans.transformComplex(Xr);
    }

//    @Override
//    public double getDX(double[] Xi, int dimTop, int dimBottom) {
//        double h = 1e-10;
//        Complex ih = new Complex(0, h);
//        Complex[] Xic = new Complex[Xi.length];
//        for (int i = 0; i < Xi.length; i++) {
//            Xic[i] = new Complex(Xi[i], 0);
//        }
//        Xic[dimBottom].add(ih);
//        Complex[] Y = getX(Xic);
//
//        return Y[dimTop].getIm() / h;
//    } 
   
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
    public double[] getXs() {
        return getX(new double[]{1.0 / 3, 1.0 / 3});
    }
}
