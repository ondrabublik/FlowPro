/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.elementType;

import flowpro.core.Parameters;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.transformation.*;

/**
 *
 * @author obublik
 */
public class triangleElementUserDef extends triangleElement {
    
    triangleElementUserDef(int order){
        super(order);
    }
    
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation2DtriangleUserDef(vertices, par);
        transform.computeTransform(vertices);
        return transform;
    }
}
