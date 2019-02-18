package flowpro.core.elementType;

import flowpro.core.Parameters;
import flowpro.core.basis.*;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.quadrature.*;
import flowpro.core.transformation.*;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class pointElement extends ElementType {
    
    public pointElement(int order) {
        this.order = order;
        nVertices = numberOfPoints();
        nFaces = numberOfEdges();
    }

    public int numberOfPoints() {
        return 1;
    }

    public int numberOfEdges() {
        return 1;
    }

    public Basis getBasis(Transformation transform) throws IOException {
        return null;
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return null;
    }

    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        return null;
    }

    public int getFaceType(int k) {
        return 1;
    }

    public int[] getFaceIndexes(int k) {
        return null;
    }
}
