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
public class lineElement extends ElementType {
    
    public lineElement(int order) {
        this.order = order;
        nVertices = numberOfPoints();
        nFaces = numberOfEdges();
    }

    public int numberOfPoints() {
        return 2;
    }

    public int numberOfEdges() {
        return 2;
    }

    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis1Dline(order);
        basis.calculateCoefficients();
        return basis;
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return qRules.qLine[order];
    }

    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv) {
        Transformation transform = new transformation1Dline();
        transform.computeTransform(vertices);
        return transform;
    }

    public int getFaceType(int k) {
        return 1;
    }

    public int[] getFaceIndexes(int k) {
        return new int[]{k};
    }
}
