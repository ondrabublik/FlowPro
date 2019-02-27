package flowpro.core.transformation;

import flowpro.api.Mat;
import flowpro.core.quadrature.Quadrature;
import java.io.Serializable;

public abstract class Transformation implements Serializable {

    public double[][] coordsX;
    public double[][] coordsXi;

    abstract public void computeTransform(double[][] Vertices);

    abstract public double[] getX(double[] Xi);

    abstract public double jacobian(double[] Xi);

    abstract public double[][] getInterpolant(Quadrature quad);

    abstract public double[] getXs();

    public double getDX(double[] Xi, int dimTop, int dimBottom) {
        double h = 1e-5;
        double[] Xih = new double[Xi.length];
        System.arraycopy(Xi, 0, Xih, 0, Xi.length);
        
        Xih[dimBottom] = Xi[dimBottom] + 4*h;
        double[] Yp4h = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] + 3*h;
        double[] Yp3h = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] + 2*h;
        double[] Yp2h = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] + h;
        double[] Yph = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] - h;
        double[] Ymh = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] - 2*h;
        double[] Ym2h = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] - 3*h;
        double[] Ym3h = getX(Xih);
        
        Xih[dimBottom] = Xi[dimBottom] - 4*h;
        double[] Ym4h = getX(Xih);
        
        double D = (-1.0/280*Yp4h[dimTop] + 4.0/105*Yp3h[dimTop] - 1.0/5*Yp2h[dimTop] + 4.0/5*Yph[dimTop] - 4.0/5*Ymh[dimTop] + 1.0/5*Ym2h[dimTop] - 4.0/105*Ym3h[dimTop] + 1.0/280*Ym4h[dimTop]) / h;
        
        return D;
    }
    
//    public double getDX(double[] Xi, int dimTop, int dimBottom) {
//        double h = 1e-5;
//        double[] Xih = new double[Xi.length];
//        System.arraycopy(Xi, 0, Xih, 0, Xi.length);
//        Xih[dimBottom] = Xi[dimBottom] + 2*h;
//        double[] Yp2h = getX(Xih);
//        
//        Xih[dimBottom] = Xi[dimBottom] + h;
//        double[] Yph = getX(Xih);
//        
//        Xih[dimBottom] = Xi[dimBottom] - h;
//        double[] Ymh = getX(Xih);
//        
//        Xih[dimBottom] = Xi[dimBottom] - 2*h;
//        double[] Ym2h = getX(Xih);
//        
//        double D = (-Yp2h[dimTop] + 8*Yph[dimTop] - 8*Ymh[dimTop] + Ym2h[dimTop])/ (12*h);
//        
//        return D;
//    }

    public double[] getXi(double[] X) {
        int dim = X.length;
        double[] Xi = new double[dim];
        double[] Xin = new double[dim];
        for (int j = 0; j < dim; j++) {
            Xi[j] = 0.5;
            Xin[j] = 0.5;
        }
        double[] V = new double[dim];
        double[][] J = new double[dim][dim];
        double h = 1e-8;
        for (int i = 0; i <= 500; i++) {
            for (int j = 0; j < dim; j++) {
                V[j] = h;
                J[j] = Mat.times(Mat.minusVec(getX(Mat.plusVec(Xi, V)), getX(Xi)), 1 / h);
                V[j] = 0;
            }
            Mat.transposeInPlace(J);
            Xin = Mat.minusVec(Xi, Mat.times(Mat.invert(J), Mat.minusVec(getX(Xi), X)));

            if (Mat.L2Norm(Mat.minusVec(Xin, Xi)) < 1e-11) {
                break;
            }
            if (i == 500) {
                System.out.println("Inverse mapping error!" + Mat.L2Norm(Mat.minusVec(Xin, Xi)));
                for (int j = 0; j < dim; j++) {
                    if (Double.isNaN(Xin[j])) {
                        Xin[j] = 0.5;
                    }
                }
                break;
            }

            System.arraycopy(Xin, 0, Xi, 0, dim);
        }
        return Xin;
    }

    public double[][] jacobiMatrix(double[] Xi) {
        int dim = Xi.length;
        double[][] J = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                J[i][j] = getDX(Xi, i, j);

            }
        }
        return J;
    }

    public double[] jacobian(Quadrature quad) {
        double[] Jac = new double[quad.nPoints];
        for (int i = 0; i < quad.nPoints; i++) {
            Jac[i] = jacobian(quad.coords[i]);
        }
        return Jac;
    }

    public double[] jacobianMatrixDirection(double[] Xi, double[] s) {
        double[][] J = jacobiMatrix(Xi);
        double[] Js = new double[Xi.length];
        for (int i = 0; i < Xi.length; i++) {
            for (int j = 0; j < Xi.length; j++) {
                Js[i] += J[i][j] * s[j];
            }
        }
        return Js;
    }

    public double[][][] transformBasis(double[][] coordsXi, double[][][] DerBasisXi) {
        double[][][] DerBasisX;
        int nPoints = DerBasisXi.length;
        int nBasis = DerBasisXi[0].length;
        int dim = DerBasisXi[0][0].length;
        DerBasisX = new double[nPoints][nBasis][dim];
        for (int p = 0; p < nPoints; p++) {
            double[][] iT = Mat.invert(jacobiMatrix(coordsXi[p]));
            for (int m = 0; m < nBasis; m++) {
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        DerBasisX[p][m][i] += iT[j][i] * DerBasisXi[p][m][j];
                    }
                }
            }
        }
        return DerBasisX;
    }

    public double[][] getX(double[][] Xi) {
        double[][] X = new double[Xi.length][Xi[0].length];
        for (int i = 0; i < Xi.length; i++) {
            X[i] = getX(Xi[i]);
        }
        return X;
    }

    public double[][] getXi(double[][] X) {
        double[][] Xi = new double[X.length][X[0].length];
        for (int i = 0; i < X.length; i++) {
            Xi[i] = getXi(X[i]);
        }
        return Xi;
    }
}
