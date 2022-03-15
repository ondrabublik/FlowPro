/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.elementType;

import flowpro.core.quadrature.Quadrature;
import flowpro.core.quadrature.QuadratureCentral;
import flowpro.core.quadrature.QuadratureLine;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author obublik
 */
public class squareElement25D extends squareElement {

    public squareElement25D(int order, int volumeQuardatureOrder, int faceQuardatureOrder) {
        super(order, volumeQuardatureOrder, faceQuardatureOrder);
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return new QuadratureSquare25D();
    }

    public class QuadratureSquare25D extends Quadrature {

        public QuadratureSquare25D() {
            dimension = 2;
            String quadratureFileName = "gaussRules/gauss/" + order + ".txt";
            QuadratureLine quad = null;
            try {
                quad = new QuadratureLine(quadratureFileName);
            } catch (IOException ex) {
                Logger.getLogger(squareElement25D.class.getName()).log(Level.SEVERE, null, ex);
            }
            int n25D = 8;
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
            System.out.println("Face weights control sum: " + controlSum);
        }

        public double[][] getCoords() {
            return coords;
        }
    }
}
