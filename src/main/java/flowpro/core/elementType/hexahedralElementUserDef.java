/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.elementType;

import flowpro.core.Parameters;
import flowpro.core.basis.*;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.transformation.*;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class hexahedralElementUserDef extends hexahedralElement {
    
    hexahedralElementUserDef(int order, int volumeQuardatureOrder, int faceQuardatureOrder){
        super(order, volumeQuardatureOrder, faceQuardatureOrder);
//        if(order == 2){
//            this.volumeQuadratureOrder = 2;
//            this.faceQuadratureOrder = 2;
//        }
    }
    
    @Override
    public Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par) {
        Transformation transform = new transformation3DhexaUserDef(par);
        transform.computeTransform(vertices);
        return transform;
    }
    
    @Override
    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis3DhexaRot(order);
        basis.calculateCoefficients();
        return basis;
    }
}
