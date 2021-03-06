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
public class hexahedralElement extends ElementType {

    public hexahedralElement(int order, int volumeQuardatureOrder, int faceQuardatureOrder) {
        this.order = order;
        this.volumeQuadratureOrder = volumeQuardatureOrder;
        this.faceQuadratureOrder = faceQuardatureOrder;
        nVertices = numberOfPoints();
        nFaces = numberOfEdges();
    }

    public int numberOfPoints() {
        return 8;
    }

    public int numberOfEdges() {
        return 6;
    }

    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis3Dhexa(order);
        basis.calculateCoefficients();
        return basis;
    }

    public Quadrature getQVolumeRule(QuadratureCentral qRules) {
        return qRules.qHexa[volumeQuadratureOrder];
    }

    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation3Dhexa();
        transform.computeTransform(vertices);
        return transform;
    }

    public int getFaceType(int k) {
        return 4;
    }

    public int[] getFaceIndexes(int k) {
        int[] faceIndexes = null;
        switch (k) {
            case 0:
                faceIndexes = new int[]{0, 1, 2, 3};
                break;
            case 1:
                faceIndexes = new int[]{4, 7, 6, 5};
                break;
            case 2:
                faceIndexes = new int[]{0, 4, 5, 1};
                break;
            case 3:
                faceIndexes = new int[]{1, 5, 6, 2};
                break;
            case 4:
                faceIndexes = new int[]{2, 6, 7, 3};
                break;
            case 5:
                faceIndexes = new int[]{3, 7, 4, 0};
                break;
        }
        return faceIndexes;
    }
}
