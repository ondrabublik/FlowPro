package flowpro.core.quadrature;

public class QuadraturePoint extends Quadrature {

    public QuadraturePoint(){
        dimension = 1;
        nPoints = 1;
        coords = new double[nPoints][dimension];
        coords[0][0] = 1;
        weights = new double[nPoints];
        weights[0] = 1;
    }

    public double[][] getCoords() {
        return coords;
    }
}
