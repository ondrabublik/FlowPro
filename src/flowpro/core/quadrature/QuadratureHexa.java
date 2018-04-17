package flowpro.core.quadrature;

import java.io.IOException;

public class QuadratureHexa extends Quadrature {

    public QuadratureHexa(String fileName) throws IOException {
        dimension = 3;
        QuadratureLine quad = new QuadratureLine(fileName);
        nPoints = quad.nPoints * quad.nPoints * quad.nPoints;
        coords = new double[nPoints][3];
        weights = new double[nPoints];

        int s = 0;
        for (int i = 0; i < quad.nPoints; i++) {
            for (int j = 0; j < quad.nPoints; j++) {
                for (int k = 0; k < quad.nPoints; k++) {
                    coords[s][0] = (1 + quad.coords[i][0]) / 2;
                    coords[s][1] = (1 + quad.coords[j][0]) / 2;
                    coords[s][2] = (1 + quad.coords[k][0]) / 2;
                    weights[s] = quad.weights[i] * quad.weights[j] * quad.weights[k] / 8;
                    ++s;
                }
            }
        }
    }

    public double[][] getCoords() {
        return coords;
    }
}
