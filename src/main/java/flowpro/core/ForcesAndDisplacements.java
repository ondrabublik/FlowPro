package flowpro.core;

import flowpro.api.FluidForces;
import flowpro.api.Mat;
import flowpro.api.MeshMove;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
public class ForcesAndDisplacements implements Serializable {

    private double dt;
    private double dto;
    private MeshMove[] mshMov;
    private FluidForces fluFor;

    public ForcesAndDisplacements(FluidForces fluFor) {
        this.fluFor = fluFor;
    }

    public ForcesAndDisplacements(double dt, double dto, MeshMove[] mshMov) {
        this.dt = dt;
        this.dto = dto;
        this.mshMov = mshMov;
    }

    public FluidForces getFluidForces() {
        return fluFor;
    }

    // sum the partial forces from subdomains
    public FluidForces combine(ForcesAndDisplacements[] forDis) {
        double[][] totalTranslationForce = null;
        double[][] totalRotationForce = null;
        double[][] boundaryForce = null; // nBoundFaces x dim
        int[][] faceIndexes = null; // element index, face index, body number
        double[][] userDef = null;

        FluidForces fluFor0 = forDis[0].getFluidForces();
        if (fluFor0.totalTranslationForce != null) {
            totalTranslationForce = Mat.copyMat(fluFor0.totalTranslationForce);
            for (int d = 1; d < forDis.length; d++) {
                FluidForces fluForAux = forDis[d].getFluidForces();
                Mat.addMat(fluForAux.totalTranslationForce, totalTranslationForce, 0);
            }
        }

        if (fluFor0.totalRotationForce != null) {
            totalRotationForce = Mat.copyMat(fluFor0.totalRotationForce);
            for (int d = 1; d < forDis.length; d++) {
                FluidForces fluForAux = forDis[d].getFluidForces();
                Mat.addMat(fluForAux.totalRotationForce, totalRotationForce, 0);
            }
        }

        if (fluFor0.userDef != null) {
            userDef = Mat.copyMat(fluFor0.userDef);
            for (int d = 1; d < forDis.length; d++) {
                FluidForces fluForAux = forDis[d].getFluidForces();
                Mat.addMat(fluForAux.userDef, userDef, 0);
            }
        }

        int s = 0;
        int dim = 0;
        for (int d = 0; d < forDis.length; d++) {
            FluidForces fluForAux = forDis[d].getFluidForces();
            if (fluForAux.boundaryForce != null) {
                s += fluForAux.boundaryForce.length;
                dim = fluForAux.boundaryForce[0].length;
            }
        }
        if (s > 0) {
            boundaryForce = new double[s][dim];
            s = 0;
            for (int d = 0; d < forDis.length; d++) {
                FluidForces fluForAux = forDis[d].getFluidForces();
                if (fluForAux.boundaryForce != null) {
                    Mat.addMat(fluForAux.boundaryForce, boundaryForce, s);
                    s += fluForAux.boundaryForce.length;
                }
            }
        }

        s = 0;
        for (int d = 0; d < forDis.length; d++) {
            FluidForces fluForAux = forDis[d].getFluidForces();
            if (fluForAux.faceIndexes != null) {
                s += fluForAux.faceIndexes.length;
            }
        }
        if (s > 0) {
            faceIndexes = new int[s][3];
            s = 0;
            for (int d = 0; d < forDis.length; d++) {
                FluidForces fluForAux = forDis[d].getFluidForces();
                if (fluForAux.faceIndexes != null) {
                    Mat.addMat(fluForAux.faceIndexes, faceIndexes, s);
                    s += fluForAux.faceIndexes.length;
                }
            }
        }
        
        return new FluidForces(totalTranslationForce, totalRotationForce, boundaryForce, faceIndexes, userDef);
    }

    public double getDt() {
        return dt;
    }

    public double getDto() {
        return dto;
    }

    public MeshMove[] getMeshMove() {
        return mshMov;
    }
}
