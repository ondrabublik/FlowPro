package flowpro.core.transformation;

import flowpro.core.quadrature.Quadrature;
import java.io.Serializable;


public abstract class FaceTransformation implements Serializable{
    
    public Transformation transform;
    
    abstract public double[] getX(double[] Xi);
    
    abstract public double[] getXiRef(double[] Xi);
    
    abstract public double jacobian(double[] Xi);
    
    abstract public double[][] getInterpolant(Quadrature quad);
    
    abstract public double[] getNormal(double[] Xi);
    
    public double[] jacobian(Quadrature quad) {
        double[] Jac = new double[quad.nPoints];
        for (int i = 0; i < quad.nPoints; i++) {
            Jac[i] = jacobian(quad.coords[i]);
        }
        return Jac;
    }
    
    public double[][] getX(double[][] Xi) {
        double[][] X = new double[Xi.length][Xi[0].length];
        for (int i = 0; i < Xi.length; i++) {
            X[i] = getX(Xi[i]);
        }
        return X;
    }
    
    public double[][] getXiRef(double[][] Xi) {
        double[][] X = new double[Xi.length][Xi[0].length];
        for (int i = 0; i < Xi.length; i++) {
            X[i] = getXiRef(Xi[i]);
        }
        return X;
    }
    
    public void setVolumeTransformation(Transformation transform){
        this.transform = transform;
    }
}
