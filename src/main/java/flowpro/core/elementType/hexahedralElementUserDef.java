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
public class hexahedralElementUserDef extends hexahedralElement {
    
    hexahedralElementUserDef(int order, int volumeQuardatureOrder, int faceQuardatureOrder){
        super(order, volumeQuardatureOrder, faceQuardatureOrder);
    }
    
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation3DhexaUserDef(par);
        transform.computeTransform(vertices);
        return transform;
    }
}
