package flowpro.core.quadrature;

import java.io.IOException;

public class QuadraturePrism extends Quadrature {

    public QuadraturePrism(Quadrature qTriangle, Quadrature qLine) throws IOException {
        dimension = 3;
        nPoints = qTriangle.nPoints * qLine.nPoints;
        coords = new double[nPoints][3];
        weights = new double[nPoints];
        
        int s = 0;
        for (int i = 0; i < qLine.nPoints; i++) {
            for (int j = 0; j < qTriangle.nPoints; j++) {
                coords[s][0] = qTriangle.coords[j][0];
                coords[s][1] = qTriangle.coords[j][1];
                coords[s][2] = (1 + qLine.coords[i][0]) / 2;
                weights[s] = qTriangle.weights[j] * qLine.weights[i] / 4;
                s++;
            }
        }
    }

    public double[][] getCoords() {
        return coords;
    }
}
