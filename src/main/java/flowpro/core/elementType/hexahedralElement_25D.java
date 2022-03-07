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
import flowpro.core.transformation.Transformation;
import flowpro.core.transformation.transformation3Dhexa;
import flowpro.core.transformation.transformation3DhexaUserDef;
import java.io.IOException;

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
        return qRules.qHexa[volumeQuadratureOrder];
    }

    @Override
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation3DhexaUserDef(par);
        transform.computeTransform(vertices);
        return transform;
    }

    @Override
    public int getFaceType(int k) {
        return 4;
    }

    @Override
    public int[] getFaceIndexes(int k) {
        int[] faceIndexes = null;
        switch (k) {
            case 0:
                faceIndexes = new int[]{0, 4, 5, 1};
                break;
            case 1:
                faceIndexes = new int[]{1, 5, 6, 2};
                break;
            case 2:
                faceIndexes = new int[]{2, 6, 7, 3};
                break;
            case 3:
                faceIndexes = new int[]{3, 7, 4, 0};
                break;
        }
        return faceIndexes;
    }
}
