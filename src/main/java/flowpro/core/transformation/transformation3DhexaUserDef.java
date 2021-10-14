package flowpro.core.transformation;

import flowpro.api.Complex;
import static flowpro.api.Complex.multiply;
import flowpro.api.DomainTransformationObject;
import flowpro.api.Mat;
import flowpro.core.Parameters;
import flowpro.core.quadrature.Quadrature;

public class transformation3DhexaUserDef extends Transformation {

    private double[][] A;
    private final double[][] invV;
    DomainTransformationObject domTrans;

    public transformation3DhexaUserDef(Parameters par) {
        domTrans = par.domainTransformationObject;
        if(domTrans == null){
            throw new UnsupportedOperationException("Domain transformation object not defined!");
        }
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

        return domTrans.transform(Xr);
    }

    public Complex[] getX(Complex[] Xi) {
        Complex[] Xr = new Complex[Xi.length];
        Complex[] Ac = new Complex[8];
        for(int i = 0; i < 8; i++){
            Ac[i] = new Complex(A[0][i], 0);
        }
        Xr[0] = new Complex();
        Xr[0].add(multiply(multiply(multiply(Xi[0],Ac[0]),Xi[1]),Xi[2]));
        Xr[0].add(multiply(multiply(Xi[0],Ac[1]),Xi[1]));
        Xr[0].add(multiply(multiply(Xi[0],Ac[2]),Xi[2]));
        Xr[0].add(multiply(multiply(Xi[1],Ac[3]),Xi[2]));
        Xr[0].add(multiply(Xi[0],Ac[4]));
        Xr[0].add(multiply(Xi[1],Ac[5]));
        Xr[0].add(multiply(Xi[2],Ac[6]));
        Xr[0].add(Ac[7]);

        for(int i = 0; i < 8; i++){
            Ac[i] = new Complex(A[1][i], 0);
        }
        Xr[1] = new Complex();
        Xr[1].add(multiply(multiply(multiply(Xi[0],Ac[0]),Xi[1]),Xi[2]));
        Xr[1].add(multiply(multiply(Xi[0],Ac[1]),Xi[1]));
        Xr[1].add(multiply(multiply(Xi[0],Ac[2]),Xi[2]));
        Xr[1].add(multiply(multiply(Xi[1],Ac[3]),Xi[2]));
        Xr[1].add(multiply(Xi[0],Ac[4]));
        Xr[1].add(multiply(Xi[1],Ac[5]));
        Xr[1].add(multiply(Xi[2],Ac[6]));
        Xr[1].add(Ac[7]);

        for(int i = 0; i < 8; i++){
            Ac[i] = new Complex(A[2][i], 0);
        }
        Xr[2] = new Complex();
        Xr[2].add(multiply(multiply(multiply(Xi[0],Ac[0]),Xi[1]),Xi[2]));
        Xr[2].add(multiply(multiply(Xi[0],Ac[1]),Xi[1]));
        Xr[2].add(multiply(multiply(Xi[0],Ac[2]),Xi[2]));
        Xr[2].add(multiply(multiply(Xi[1],Ac[3]),Xi[2]));
        Xr[2].add(multiply(Xi[0],Ac[4]));
        Xr[2].add(multiply(Xi[1],Ac[5]));
        Xr[2].add(multiply(Xi[2],Ac[6]));
        Xr[2].add(Ac[7]);

        return domTrans.transformComplex(Xr);
    }

//    @Override
//    public double getDX(double[] Xi, int dimTop, int dimBottom) {
//        double h = 1e-8;
//        Complex ih = new Complex(0, h);
//        Complex[] Xic = new Complex[Xi.length];
//        for (int i = 0; i < Xi.length; i++) {
//            Xic[i] = new Complex(Xi[i], 0);
//        }
//        Xic[dimBottom].add(ih);
//        Complex[] Y = getX(Xic);
//        return Y[dimTop].getIm() / h;
//    }

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
