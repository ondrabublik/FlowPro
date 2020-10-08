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

    public double[][] totalTranslationForce;
    public double[][] totalRotationForce;
    public double[][] userDef;

    public Deformation2DRigid(Parameters par, Equation eqn, int[][] TEale) throws IOException {
        super(par, eqn, TEale);
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
        for (int i = 0; i < elems.length; i++) {
            for (int j = 0; j < elems[i].nVertices; j++) {
                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].vertices[j][d] = 0;
                }
                for (int k = 0; k < nBodies; k++) {
                    double[] moveTranslation = mshMov[k].getTranslation();
                    double[] moveRotation = mshMov[k].getRotation();
                    double xNew = (Math.cos(moveRotation[0]) * (elems[i].vertices0[j][0] - center[0][k]) - Math.sin(moveRotation[0]) * (elems[i].vertices0[j][1] - center[1][k]) + center[0][k] - elems[i].vertices0[j][0]) * nBodies + elems[i].vertices0[j][0] + moveTranslation[0] * nBodies;
                    double yNew = (Math.sin(moveRotation[0]) * (elems[i].vertices0[j][0] - center[0][k]) + Math.cos(moveRotation[0]) * (elems[i].vertices0[j][1] - center[1][k]) + center[1][k] - elems[i].vertices0[j][1]) * nBodies + elems[i].vertices0[j][1] + moveTranslation[1] * nBodies;
                    elems[i].vertices[j][0] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][0] + elems[i].blendFun[j][k] * xNew) / nBodies;
                    elems[i].vertices[j][1] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][1] + elems[i].blendFun[j][k] * yNew) / nBodies;
                }

                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].U[j][d] = a1 * elems[i].vertices[j][d] + a2 * elems[i].verticesOld[j][d] + a3 * elems[i].verticesOld2[j][d];
                }
            }
        }
    }

	@Override
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
        totalTranslationForce = new double[2][nBodies];
        totalRotationForce = new double[1][nBodies];
        for (int b = 0; b < nBodies; b++) {
            double[] moveTranslation = mshMov[b].getTranslation();
            for (Element elem : elems) {
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TEale[k] == b + 2 && elem.insideMetisDomain) {
						Face face = elem.Int.faces[k];
                        double[] Jac = face.JacobianFace;
                        double[] weights = face.weightsFace;
//                        double[][] baseLeft = face.basisFaceLeft;
//						double[][][] derBaseLeft = face.dXbasisFaceLeft;
                        double fx = 0;
                        double fy = 0;
                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
							double[] wL = face.evalWLeft(elem.W, p);
							double[] derWL = face.evalDerWLeft(elem.W, p);
//                            double[] WL = new double[elem.getNEqs()];
//                            for (int j = 0; j < elem.getNEqs(); j++) {
//                                for (int m = 0; m < elem.nBasis; m++) {
//                                    WL[j] = WL[j] + elem.W[j * elem.nBasis + m] * baseLeft[p][m];
//                                }
//                            }
//                            double pressure = eqn.pressure(WL);
							double[] normalStress = eqn.normalStress(wL, derWL, elem.n[k][p]);
//                            fx += Jac[p] * weights[p] * elem.n[k][p][0] * pressure;
//                            fy += Jac[p] * weights[p] * elem.n[k][p][1] * pressure;
							fx -= Jac[p] * weights[p] * normalStress[0];
                            fy -= Jac[p] * weights[p] * normalStress[1];
                        }
                        totalTranslationForce[0][b] += fx;
                        totalTranslationForce[1][b] += fy;
                        totalRotationForce[0][b] += -fx * (elem.Xes[k][1] - (center[1][b] + moveTranslation[1])) + fy * (elem.Xes[k][0] - (center[0][b] + moveTranslation[0]));
                    }
                }
            }
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
    public FluidForces getFluidForces() {
        return new FluidForces(totalTranslationForce, totalRotationForce, null, null, userDef);
    }
}
