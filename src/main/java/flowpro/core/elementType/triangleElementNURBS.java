/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.elementType;

import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.transformation.*;

/**
 *
 * @author obublik
 */
public class triangleElementNURBS extends triangleElement {
    
    triangleElementNURBS(int order){
        super(order);
    }
    
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv) {
        Transformation transform = new transformation2DtriangleNURBS(vertices);
        transform.computeTransform(vertices);
        return transform;
    }
}
