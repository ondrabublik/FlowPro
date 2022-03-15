package flowpro.core.quadrature;

import java.io.IOException;

public class QuadratureSquare25D extends Quadrature {

    public QuadratureSquare25D(int order) throws IOException {
        dimension = 2;
        
        String quadratureFileName = "gaussRules/gauss/" + order + ".txt";
        QuadratureLine quad = null;
        quad = new QuadratureLine(quadratureFileName);

        int n25D = 12;
        nPoints = quad.nPoints * n25D;
        coords = new double[nPoints][2];
        weights = new double[nPoints];

        double controlSum = 0;
        int s = 0;
        double dx = 1.0 / n25D;
        for (int i = 0; i < quad.nPoints; i++) {
            for (int j = 0; j < n25D; j++) {
                coords[s][0] = (1 + quad.coords[i][0]) / 2;
                coords[s][1] = dx / 2 + dx * j;
                weights[s] = quad.weights[i] / 2 * 1.0 / n25D;
                controlSum += weights[s];
                ++s;
            }
        }
        //System.out.println(controlSum);
    }

    public double[][] getCoords() {
        return coords;
    }
}
