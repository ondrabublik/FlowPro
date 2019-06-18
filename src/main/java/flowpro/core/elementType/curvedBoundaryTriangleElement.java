package flowpro.core.elementType;

import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.transformation.*;

/**
 *
 * @author obublik
 */
public class curvedBoundaryTriangleElement extends triangleElement{
    
    curvedBoundaryTriangleElement(int order, int volumeQuardatureOrder, int faceQuardatureOrder){
            super(order, volumeQuardatureOrder, faceQuardatureOrder);
    }
    
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv) {
        Transformation transform = new transformation2DTriangleCurved(fCurv);
        transform.computeTransform(vertices);
        return transform;
    }
}
