/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.elementType;

import flowpro.core.Parameters;
import flowpro.core.basis.Basis;
import flowpro.core.basis.basis3D_25D;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.quadrature.Quadrature;
import flowpro.core.quadrature.QuadratureCentral;
import flowpro.core.quadrature.QuadratureLine;
import flowpro.core.transformation.Transformation;
import flowpro.core.transformation.transformation3DhexaUserDef;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author obublik
 */
public class hexahedralElement_25D extends ElementType {

    public hexahedralElement_25D(int order, int volumeQuardatureOrder, int faceQuardatureOrder) {
        this.order = order;
        this.volumeQuadratureOrder = volumeQuardatureOrder;
        this.faceQuadratureOrder = faceQuardatureOrder;
        nVertices = numberOfPoints();
        nFaces = numberOfEdges();
    }

    @Override
    public int numberOfPoints() {
        return 8;
    }

    @Override
    public int numberOfEdges() {
        return 4;
    }

    @Override
    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis3D_25D(order);
        basis.calculateCoefficients();
        return basis;
    }

    @Override
    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return new Quadrature25D();
    }

    @Override
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation3DhexaUserDef(par);
        transform.computeTransform(vertices);
        return transform;
    }

    @Override
    public int getFaceType(int k) {
        return 43;
    }

    @Override
    public int[] getFaceIndexes(int k) {
        int[] faceIndexes = null;
        switch (k) {
            case 0:
                faceIndexes = new int[]{1, 0, 4, 5};
                break;
            case 1:
                faceIndexes = new int[]{2, 1, 5, 6};
                break;
            case 2:
                faceIndexes = new int[]{3, 2, 6, 7};
                break;
            case 3:
                faceIndexes = new int[]{0, 3, 7, 4};
                break;
        }
        return faceIndexes;
    }

    public class Quadrature25D extends Quadrature {

        
        public Quadrature25D() {
            dimension = 3;
            String quadratureFileName = "gaussRules/gauss/" + volumeQuadratureOrder + ".txt";
            QuadratureLine quad = null;
            try {
                quad = new QuadratureLine(quadratureFileName);
            } catch (IOException ex) {
                Logger.getLogger(hexahedralElement_25D.class.getName()).log(Level.SEVERE, null, ex);
            }
            int n25D = 12;
            nPoints = quad.nPoints * quad.nPoints * n25D;
            coords = new double[nPoints][3];
            weights = new double[nPoints];

            double controlSum = 0;
            int s = 0;
            double dx = 1.0 / n25D;
            for (int i = 0; i < quad.nPoints; i++) {
                for (int j = 0; j < quad.nPoints; j++) {
                    for (int k = 0; k < n25D; k++) {
                        coords[s][0] = (1 + quad.coords[i][0]) / 2;
                        coords[s][1] = (1 + quad.coords[j][0]) / 2;
                        coords[s][2] = dx / 2 + dx * k;
                        weights[s] = quad.weights[i] * quad.weights[j] / 4 * 1.0 / n25D;
                        controlSum += weights[s];
                        ++s;
                    }
                }
            }
            //System.out.println("Volume weights control sum: " + controlSum);
        }

        public double[][] getCoords() {
            return coords;
        }
    }
}
