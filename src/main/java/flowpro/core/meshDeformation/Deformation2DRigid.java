package flowpro.core.meshDeformation;

import flowpro.api.Equation;
import flowpro.api.FluidForces;
import flowpro.api.MeshMove;
import flowpro.core.Integration.Face;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import java.io.*;

/**
 *
 * @author obublik
 */
public class Deformation2DRigid extends Deformation {

    private final FluidForces[] fluidForces;

    public Deformation2DRigid(Parameters par, Equation eqn, int[][] TEale) throws IOException {
        super(par, eqn, TEale);

        fluidForces = new FluidForces[nBodies];
    }

    @Override
    public void newMeshPositionAndVelocity(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {

        double a1 = 1.0 / dt;
        double a2 = -1.0 / dt;
        double a3 = 0;

        if (timeOrder > 1) {
            a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
            a2 = -(dt + dto) / (dt * dto);  // -2/dt;
            a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
        }

        // multiple blending function
        for (Element elem : elems) {
            for (int k = 0; k < nBodies; k++) {
                double[] moveTranslation = mshMov[k].getTranslation();
                double[] moveRotation = mshMov[k].getRotation();
                for (int j = 0; j < elem.nVertices; j++) {
                    for (int d = 0; d < elem.dim; d++) {
                        elem.vertices[j][d] = 0;
                    }
                    double xNew = (Math.cos(moveRotation[0]) * (elem.vertices0[j][0] - center[0][k]) - Math.sin(moveRotation[0]) * (elem.vertices0[j][1] - center[1][k]) + center[0][k] - elem.vertices0[j][0]) * nBodies + elem.vertices0[j][0] + moveTranslation[0] * nBodies;
                    double yNew = (Math.sin(moveRotation[0]) * (elem.vertices0[j][0] - center[0][k]) + Math.cos(moveRotation[0]) * (elem.vertices0[j][1] - center[1][k]) + center[1][k] - elem.vertices0[j][1]) * nBodies + elem.vertices0[j][1] + moveTranslation[1] * nBodies;
                    elem.vertices[j][0] += ((1 - elem.blendFun[j][k]) * elem.vertices0[j][0] + elem.blendFun[j][k] * xNew) / nBodies;
                    elem.vertices[j][1] += ((1 - elem.blendFun[j][k]) * elem.vertices0[j][1] + elem.blendFun[j][k] * yNew) / nBodies;
                }
            }
        }

        for (Element elem : elems) {
            for (int j = 0; j < elem.nVertices; j++) {
                for (int d = 0; d < elem.dim; d++) {
                    elem.U[j][d] = a1 * elem.vertices[j][d] + a2 * elem.verticesOld[j][d] + a3 * elem.verticesOld2[j][d];
                }
            }
        }
    }

    @Override
    public void nextTimeLevel(Element[] elems) {
        for (Element elem : elems) {
            for (int j = 0; j < elem.nVertices; j++) {
                for (int d = 0; d < elem.dim; d++) {
                    elem.verticesOld2[j][d] = elem.verticesOld[j][d];
                    elem.verticesOld[j][d] = elem.vertices[j][d];
                }
            }
        }
    }

//	private double[] evalW(Element elem, double[] base) {
//		double[] W = new double[elem.getNEqs()];
//		for (int m = 0; m < elem.getNEqs(); m++) {
//			for (int i = 0; i < elem.nBasis; i++) {
//				W[m] = W[m] + elem.W[m * elem.nBasis + i] * base[i];
//			}
//		}		
//		return W;
//	}
//	
//	private double[] derEvalW(Element elem, double[][] derBase) {
//		int dim = derBase[0].length;
//		int nEqs = elem.getNEqs();
//		double[] derW = new double[nEqs * dim];
//		for (int m = 0; m < elem.getNEqs(); m++) {
//			for (int i = 0; i < elem.nBasis; i++) {
//				for (int d = 0; d < dim; ++d) {
//					derW[nEqs * d + m] += elem.W[m * elem.nBasis + i] * derBase[i][d];
//				}
//			}
//		}		
//		return derW;
//	}
    @Override
    public void calculateForces(Element[] elems, MeshMove[] mshMov) {

        for (int b = 0; b < nBodies; b++) {
            double[] moveTranslation = mshMov[b].getTranslation();
            double[] totalTranslationForce = new double[2];
            double[] totalRotationForce = new double[1];
            for (Element elem : elems) {
                for (int r = 0; r < elem.nFaces; r++) {
                    if (elem.TEale[r] == b + 2 && elem.insideMetisDomain) {
                        Face face = elem.Int.faces[r];

                        double[] force = face.integrateLeft((double[] w, double[] dw, double[] n)
                                -> eqn.stressVector(w, dw, n), elem.W, elem.n[r]);
                        double fx = -force[0];
                        double fy = -force[1];

//						double[] Jac = face.JacobianFace;
//                        double[] weights = face.weightsFace;
//                        double fx = 0;
//                        double fy = 0;
//                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
//							double[] wL = face.evalWLeft(elem.W, p);
//							double[] derWL = face.evalDerWLeft(elem.W, p);
//							double[] normalStress = eqn.normalStress(wL, derWL, elem.n[k][p]);
//							fx -= Jac[p] * weights[p] * normalStress[0];
//                            fy -= Jac[p] * weights[p] * normalStress[1];
//                        }
                        totalTranslationForce[0] += fx;
                        totalTranslationForce[1] += fy;
                        totalRotationForce[0] += -fx * (elem.Xes[r][1] - (center[1][b] + moveTranslation[1])) + fy * (elem.Xes[r][0] - (center[0][b] + moveTranslation[0]));
                    }
                }
            }
            fluidForces[b] = new FluidForces(totalTranslationForce, totalRotationForce);
        }

        // user defined totalTranslationForce term
//        userDef = new double[nUserValues][nBodies];
//        for (int b = 0; b < nBodies; b++) {
//            for (Element elem : elems) {
//                for (int k = 0; k < elem.nFaces; k++) {
//                    if (elem.TEale[k] == b + 2 && elem.insideMetisDomain) {
//                        double[] Jac = elem.Int.faces[k].JacobianFace;
//                        double[] weights = elem.Int.faces[k].weightsFace;
//                        double[][] baseLeft = elem.Int.faces[k].basisFaceLeft;
//                        double fyTop = 0;
//                        double fyBottom = 0;
//                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
////                            double[] WL = new double[elem.getNEqs()];
////                            for (int j = 0; j < elem.getNEqs(); j++) {
////                                for (int m = 0; m < elem.nBasis; m++) {
////                                    WL[j] = WL[j] + elem.W[j * elem.nBasis + m] * baseLeft[p][m];
////                                }
////                            }
////                            double pressure = eqn.pressure(WL);
//							double[] wL = face.evalWLeft(elem.W, p);
//							double[] derWL = face.evalDerWLeft(elem.W, p);
//							double[] normalStress = eqn.normalStress(WL, WL, elem.n[k][p]);
//                            if (elem.n[k][p][1] > 0) {
////                                fyTop += Jac[p] * weights[p] * elem.n[k][p][1] * pressure;
//								fyTop -= Jac[p] * weights[p] * normalStress[1];
//                            } else {
////                                fyBottom += Jac[p] * weights[p] * elem.n[k][p][1] * pressure;
//								fyBottom -= Jac[p] * weights[p] * normalStress[1];
//                            }
//                        }
//                        userDef[0][b] += fyTop;
//                        userDef[1][b] += fyBottom;
//                    }
//                }
//            }
//        }
    }

    @Override
    public FluidForces[] getFluidForces() {
        return fluidForces;
    }
}
