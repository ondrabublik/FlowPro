package flowpro.core.elementType;

import flowpro.core.basis.*;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.quadrature.*;
import flowpro.core.transformation.*;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class prismElement extends ElementType {

    public prismElement(int order) {
        this.order = order;
        nVertices = numberOfPoints();
        nFaces = numberOfEdges();
    }

    public int numberOfPoints() {
        return 6;
    }

    public int numberOfEdges() {
        return 5;
    }

    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis3Dprism(order);
        basis.calculateCoefficients();
        return basis;
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return qRules.qPrism[order];
    }

    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv) {
        Transformation transform = new transformation3Dprism();
        transform.computeTransform(vertices);
        return transform ;
    }

    public int getFaceType(int k) {
        if (k < 2) {
            return 3;
        } else {
            return 4;
        }
    }

    public int[] getFaceIndexes(int k) {
        int[] faceIndexes = null;
        switch (k) {
            case 0:
                faceIndexes = new int[]{0, 1, 2};
                break;
            case 1:
                faceIndexes = new int[]{3, 5, 4};
                break;
            case 2:
                faceIndexes = new int[]{0, 3, 4, 1};
                break;
            case 3:
                faceIndexes = new int[]{1, 4, 5, 2};
                break;
            case 4:
                faceIndexes = new int[]{2, 5, 3, 0};
                break;
        }
        return faceIndexes;
    }
}
