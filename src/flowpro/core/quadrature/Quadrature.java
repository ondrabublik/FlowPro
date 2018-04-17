package flowpro.core.quadrature;

import java.io.Serializable;


/**
 * Quadrature rules.
 *
 * @author ales
 */
public abstract class Quadrature implements Serializable {
    
    public int dimension;     // dimension
    public int nPoints;       // number of quadrature coords
    public double[][] coords; // quadrature coords
    public double[] weights;  // quadrature weights
    
    abstract public double[][] getCoords();
}