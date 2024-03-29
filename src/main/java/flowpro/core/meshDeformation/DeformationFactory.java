package flowpro.core.meshDeformation;

import flowpro.api.Equation;
import flowpro.core.Parameters;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class DeformationFactory {

    public Deformation getDeformation(Parameters par, Equation eqn, double[][] PXY, int[][] TP, int[][] TEale) throws IOException {
        switch (par.meshDeformationType) {
            case "Rigid2D":
                return new Deformation2DRigid(par, eqn, TEale);
            case "Elastic2D":
                return new Deformation2DElastic1(par, eqn, PXY, TP, TEale);
            case "Rigid3D":
                return new Deformation3DRigid(par, eqn, TEale);
            case "Elastic3D":
                return null;
            case "Elastic3DBeam":
                return new Deformation3DBeam(par, eqn, TEale);
            case "Test2D":
                return new Deformation2DTest(par, eqn, TEale);
            default:
                System.out.println("Mesh deformation algorithm not found:" + par.meshDeformationType);
                return null;
        }

    }
}
