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
public class tetrahedralElement extends ElementType {

    public tetrahedralElement(int order, int volumeQuardatureOrder, int faceQuardatureOrder) {
        this.order = order;
        this.volumeQuadratureOrder = volumeQuardatureOrder;
        this.faceQuadratureOrder = faceQuardatureOrder;
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
        Basis basis = new basis3Dtetra(order);
        basis.calculateCoefficients();
        return basis;
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return qRules.qTetra[volumeQuadratureOrder];
    }

    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation3Dtetra(vertices);
        transform.computeTransform(vertices);
        return transform;
    }

    public int getFaceType(int k) {
        return 3;
    }

    public int[] getFaceIndexes(int k) {
        int[] faceIndexes = null;
        switch (k) {
            case 0:
                faceIndexes = new int[]{0, 2, 1};
                break;
            case 1:
                faceIndexes = new int[]{0, 1, 3};
                break;
            case 2:
                faceIndexes = new int[]{1, 2, 3};
                break;
            case 3:
                faceIndexes = new int[]{2, 0, 3};
                break;
        }
        return faceIndexes;
    }
}
