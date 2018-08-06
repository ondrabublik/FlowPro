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
public class squareElement extends ElementType{
    
    public squareElement(int order) {
        this.order = order;
        nVertices = numberOfPoints();
        nFaces = numberOfEdges();
    }

    public int numberOfPoints() {
        return 4;
    }

    public int numberOfEdges() {
        return 4;
    }

    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis2Dsquare(order);
        basis.calculateCoefficients();
        return basis;
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return qRules.qSquare[order];
    }

    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv) {
        Transformation transform = new transformation2Dsquare();
        transform.computeTransform(vertices);
        return transform;
    }

    public int getFaceType(int k) {
        return 2;
    }

    public int[] getFaceIndexes(int k) {
        return new int[]{k,(k + 1) % nVertices};
    }
}
