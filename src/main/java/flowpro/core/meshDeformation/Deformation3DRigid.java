package flowpro.core.meshDeformation;

import flowpro.api.Equation;
import flowpro.api.FluidForces;
import flowpro.api.Mat;
import flowpro.api.MeshMove;
import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import java.io.*;

/**
 *
 * @author obublik
 */
public class Deformation3DRigid extends Deformation {

    public double[][] totalTranslationForce;
    public double[][] totalRotationForce;

    public Deformation3DRigid(Parameters par, Equation eqn, int[][] TEale) throws IOException {
        super(par, eqn, TEale);
    }

    public void newMeshPosition(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {

        double a1 = 1.0 / dt;
        double a2 = -1.0 / dt;
        double a3 = 0;

        if (timeOrder > 1) {
            a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
            a2 = -(dt + dto) / (dt * dto);  // -2/dt;
            a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
        }

        // multiple blending function
        for (int i = 0; i < elems.length; i++) {
            for (int j = 0; j < elems[i].nVertices; j++) {
                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].vertices[j][d] = 0;
                }
                for (int k = 0; k < nBodies; k++) {
                    double[] moveTranslation = mshMov[k].getTranslation();
                    double[] moveRotation = mshMov[k].getRotation();
                    double[][] Rotmat = rotationMatrix(moveRotation);
                    double[] Xcenter = new double[]{center[0][k], center[1][k], center[2][k]};
                    double[] rotation = Mat.times(Rotmat, Mat.minusVec(elems[i].vertices0[j], Xcenter));
                    double xNew = (rotation[0] + Xcenter[0] - elems[i].vertices0[j][0]) * nBodies + elems[i].vertices0[j][0] + moveTranslation[0] * nBodies;
                    double yNew = (rotation[1] + Xcenter[1] - elems[i].vertices0[j][1]) * nBodies + elems[i].vertices0[j][1] + moveTranslation[1] * nBodies;
                    double zNew = (rotation[2] + Xcenter[2] - elems[i].vertices0[j][2]) * nBodies + elems[i].vertices0[j][2] + moveTranslation[2] * nBodies;
                    elems[i].vertices[j][0] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][0] + elems[i].blendFun[j][k] * xNew) / nBodies;
                    elems[i].vertices[j][1] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][1] + elems[i].blendFun[j][k] * yNew) / nBodies;
                    elems[i].vertices[j][2] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][2] + elems[i].blendFun[j][k] * zNew) / nBodies;
                }

                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].U[j][d] = a1 * elems[i].vertices[j][d] + a2 * elems[i].verticesOld[j][d] + a3 * elems[i].verticesOld2[j][d];
                }
            }
        }
    }

    public void nextTimeLevel(Element[] elems) {
        for (int i = 0; i < elems.length; i++) {
            for (int j = 0; j < elems[i].nVertices; j++) {
                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].verticesOld2[j][d] = elems[i].verticesOld[j][d];
                    elems[i].verticesOld[j][d] = elems[i].vertices[j][d];
                }
            }
        }
    }

    public void calculateForces(Element[] elems, MeshMove[] mshMov) {
        totalTranslationForce = new double[3][nBodies];
        totalRotationForce = new double[3][nBodies];
        for (int b = 0; b < nBodies; b++) {
            double[] moveTranslation = mshMov[b].getTranslation();
            for (Element elem : elems) {
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TEale[k] == b + 2 && elem.insideMetisDomain) {
                        double[] Jac = elem.Int.faces[k].JacobianFace;
                        double[] weights = elem.Int.faces[k].weightsFace;
                        double[][] baseLeft = elem.Int.faces[k].basisFaceLeft;
                        double[] f = new double[3];
                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
                            double[] WL = new double[elem.getNEqs()];
                            for (int j = 0; j < elem.getNEqs(); j++) {
                                for (int m = 0; m < elem.nBasis; m++) {
                                    WL[j] = WL[j] + elem.W[j * elem.nBasis + m] * baseLeft[p][m];
                                }
                            }
                            double pressure = eqn.pressure(WL);
                            for (int d = 0; d < elem.dim; d++) {
                                f[d] += Jac[p] * weights[p] * elem.n[k][p][d] * pressure;
                            }
                        }
                        double[] Xcenter = new double[]{center[0][b], center[1][b], center[2][b]};
                        double[] M = Mat.cross(f, Mat.minusVec(elem.Xes[k], Mat.plusVec(Xcenter,moveTranslation)));
                        for (int d = 0; d < elem.dim; d++) {
                            totalTranslationForce[d][b] += f[d];
                            totalRotationForce[d][b] += M[d];
                        }
                    }
                }
            }
        }
    }

    public FluidForces getFluidForces() {
        return new FluidForces(totalTranslationForce, totalRotationForce, null, null, null);
    }

    double[][] rotationMatrix(double[] angles) {
        double[][] Ra = new double[][]{{1, 0, 0}, {0, Math.cos(angles[0]), -Math.sin(angles[0])}, {0, Math.sin(angles[0]), Math.cos(angles[0])}};
        double[][] Rb = new double[][]{{Math.cos(angles[1]), 0, -Math.sin(angles[1])}, {0, 1, 0}, {Math.sin(angles[1]), 0, Math.cos(angles[1])}};
        double[][] Rc = new double[][]{{Math.cos(angles[2]), -Math.sin(angles[2]), 0}, {Math.sin(angles[2]), Math.cos(angles[2]), 0}, {0, 0, 1}};
        
        return Mat.times(Ra, Mat.times(Rb, Rc));
    }
}
