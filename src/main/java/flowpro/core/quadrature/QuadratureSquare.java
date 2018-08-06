package flowpro.core.quadrature;

import java.io.IOException;

public class QuadratureSquare extends Quadrature {
    
    public QuadratureSquare(String fileName) throws IOException {
        dimension = 2;
        QuadratureLine quad = new QuadratureLine(fileName);
        nPoints = quad.nPoints * quad.nPoints;
        coords = new double[nPoints][2];
        weights = new double[nPoints];

        int s = 0;
        for (int i = 0; i < quad.nPoints; i++) {
            for (int j = 0; j < quad.nPoints; j++) {
                coords[s][0] = (1 + quad.coords[i][0]) / 2;
                coords[s][1] = (1 + quad.coords[j][0]) / 2;
                weights[s] = quad.weights[i] * quad.weights[j] / 4;
                ++s;
            }
        }
    }

    public double[][] getCoords() {
        return coords;
    }
}
