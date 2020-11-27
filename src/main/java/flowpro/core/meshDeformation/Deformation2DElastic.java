package flowpro.core.meshDeformation;

import flowpro.api.Equation;
import flowpro.api.FluidForces;
import flowpro.api.MeshMove;
import flowpro.core.Integration;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class Deformation2DElastic extends Deformation {

	public double[][] boundaryForce;  // smazat
	public int[][] faceIndexes;  // smazat
	
	private final FluidForces[] fluidForces;

	public Deformation2DElastic(Parameters par, Equation eqn, int[][] TEale) throws IOException {
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
		for (int i = 0; i < elems.length; i++) {
			for (int j = 0; j < elems[i].nVertices; j++) {
				for (int d = 0; d < elems[i].dim; d++) {
					elems[i].vertices[j][d] = 0;
				}
				for (int k = 0; k < nBodies; k++) {
					double[] moveTranslation = mshMov[k].getTranslation();
					double[] moveRotation = mshMov[k].getRotation();
					double[][] rbfCoeffs = mshMov[k].getDeformationCoefficients();
					double[] Xdef = new double[elems[i].dim];
					if (rbfCoeffs != null) {
						double[][] boundaryPointsCoords = mshMov[k].getboundaryPointsCoords();
						for (int p = 0; p < rbfCoeffs[0].length; p++) {
							double rbf = radialBasisFunction(elems[i].vertices0[j], boundaryPointsCoords[p]);
							for (int d = 0; d < elems[i].dim; d++) {
								Xdef[d] += rbfCoeffs[d][p] * rbf;
							}
						}
					}
					double xNew = (Math.cos(moveRotation[0]) * (elems[i].vertices0[j][0] - center[0][k]) - Math.sin(moveRotation[0]) * (elems[i].vertices0[j][1] - center[1][k]) + center[0][k] - elems[i].vertices0[j][0]) * nBodies + elems[i].vertices0[j][0] + moveTranslation[0] * nBodies + Xdef[0] * nBodies;
					double yNew = (Math.sin(moveRotation[0]) * (elems[i].vertices0[j][0] - center[0][k]) + Math.cos(moveRotation[0]) * (elems[i].vertices0[j][1] - center[1][k]) + center[1][k] - elems[i].vertices0[j][1]) * nBodies + elems[i].vertices0[j][1] + moveTranslation[1] * nBodies + Xdef[1] * nBodies;
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

//	@Override
	public void calculateForces1(Element[] elems, MeshMove[] mshMov) {
		boundaryForce = null;
		faceIndexes = null;

		int s = 0;
		for (Element elem : elems) {
			for (int k = 0; k < elem.nFaces; k++) {
				if (elem.TEale[k] > 1 && elem.insideMetisDomain) {
					s += elem.Int.faces[k].nIntEdge;
				}
			}
		}

		if (s > 0) {
			int dim = elems[0].dim;
			boundaryForce = new double[s][2 * dim]; // s x [Fx,Fy,Fz,xes,yes,zes]
			faceIndexes = new int[s][3];
			s = 0;
			for (int k = 0; k < elems.length; k++) {
				Element elem = elems[k];
				for (int r = 0; r < elem.nFaces; r++) {
					if (elem.TEale[r] > 1 && elem.insideMetisDomain) {
						Integration.Face face = elem.Int.faces[r];
						for (int p = 0; p < face.nIntEdge; p++) { // edge integral                            
							double[] WL = face.evalWLeft(elem.W, p);
							double pressure = eqn.pressure(WL);
//							double[] derWL = face.evalDerWLeft(elem.W, p);
//							double[] stressVector = eqn.stressVector(WL, derWL, elem.n[r][p]);
							for (int d = 0; d < dim; d++) {
                                boundaryForce[s][d] = elem.n[r][p][d] * pressure;
//								boundaryForce[s][d] = -stressVector[d];
							}

							System.arraycopy(elem.Xes[r], 0, boundaryForce[s], dim, dim);
							faceIndexes[s][0] = k;
							faceIndexes[s][1] = r;
							faceIndexes[s][2] = elem.TEale[r];

							s++;
						}
					}
				}
			}
		}
	}
		
//	@Override
	public void calculateForcesOld(Element[] elems, MeshMove[] mshMov) {
		boundaryForce = null;
		faceIndexes = null;

		int s = 0;
		for (Element elem : elems) {
			for (int k = 0; k < elem.nFaces; k++) {
				if (elem.TEale[k] > 1 && elem.insideMetisDomain) {
					s++;
				}
			}
		}
		/*
        if (s > 0) {
            int dim = elems[0].dim;
            boundaryForce = new double[s][2*dim]; // s x [Fx,Fy,Fz,xes,yes,zes]
            faceIndexes = new int[s][3];
            s = 0;
            for (int i = 0; i < elems.length; i++) {
                Element elem = elems[i];
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TEale[k] > 1 && elem.insideMetisDomain) {
                        double[] Jac = elem.Int.faces[k].JacobianFace;
                        double[] weights = elem.Int.faces[k].weightsFace;
                        double[][] baseLeft = elem.Int.faces[k].basisFaceLeft;
                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
                            double[] WL = new double[elem.getNEqs()];
                            for (int j = 0; j < elem.getNEqs(); j++) {
                                for (int m = 0; m < elem.nBasis; m++) {
                                    WL[j] = WL[j] + elem.W[j * elem.nBasis + m] * baseLeft[p][m];
                                }
                            }
                            double pressure = eqn.pressure(WL);
                            for (int d = 0; d < dim; d++) {
                                boundaryForce[s][d] += Jac[p] * weights[p] * elem.n[k][p][d] * pressure;
                            }
                        }
                        System.arraycopy(elem.Xes[k], 0, boundaryForce[s], dim, dim);
                        faceIndexes[s][0] = i;
                        faceIndexes[s][1] = k;
                        faceIndexes[s][2] = elem.TEale[k];
                        s++;
                    }
                }
            }
        }*/

		if (s > 0) {
			int dim = elems[0].dim;
			// vektor napeti a souradnice stredu strany
			boundaryForce = new double[s][2 * dim]; // [sigmaX, sigmaY, sigmaZ, xes, yes, zes]
			faceIndexes = new int[s][3];
			s = 0;
			for (int k = 0; k < elems.length; k++) {
				Element elem = elems[k];
				for (int r = 0; r < elem.nFaces; r++) {
					if (elem.TEale[r] > 1 && elem.insideMetisDomain) {
//						Integration.Face face = elem.Int.faces[r];
//						double[] WL = face.evalWLeft(elem.W, 0);
//						double[] derWL = face.evalDerWLeft(elem.W, 0);
//						double[] stressVector = eqn.stressVector(WL, derWL, elem.n[r][0]);
						double pressure = eqn.pressure(elem.calculateAvgW());
						for (int d = 0; d < dim; d++) {
//							boundaryForce[s][d] = pressure;
                            boundaryForce[s][d] = pressure * elem.n[r][0][d];  // prumernou normalu nebo integrovat!!!
//							boundaryForce[s][d] = -stressVector[d];
						}
						System.arraycopy(elem.Xes[r], 0, boundaryForce[s], dim, dim);
						faceIndexes[s][0] = k;
						faceIndexes[s][1] = r;
						faceIndexes[s][2] = elem.TEale[r];
						s++;
					}
				}
			}
		}
	}

	@Override
	public FluidForces[] getFluidForces() {
		return fluidForces;
	}

	public double radialBasisFunction(double[] a, double[] b) {
		double norm = 0;
		for (int i = 0; i < a.length; i++) {
			norm += (a[i] - b[i]) * (a[i] - b[i]);
		}
		norm = Math.sqrt(norm);
		return 1 - norm;
	}
	
	@Override
	public void calculateForces(Element[] elems, MeshMove[] mshMov) {
		BoundaryIndices[] boundaryIndexesOfBodies = getBoundaryIndexesOfBodies(elems);
		
		int dim = elems[0].dim;
		for (int b = 0; b < nBodies; b++) {	
			BoundaryIndices boundaryIndexes = boundaryIndexesOfBodies[b];
			
			double[] force = new double[dim];
			double[][] stressVectors = new double[boundaryIndexes.nIntegrationPoints][dim];
			double[][] stressVectorPositions = new double[boundaryIndexes.nIntegrationPoints][dim];
			
			for (int i = 0; i < boundaryIndexes.elemInds.length; i++) {
				Element elem = elems[boundaryIndexes.elemInds[i]];
				int faceIdx = boundaryIndexes.faceInds[i];
				Integration.Face face = elem.Int.faces[faceIdx];
				for (int p = 0; p < face.nIntEdge; p++) {
					int globalIdx = i * face.nIntEdge + p;
					double[] WL = face.evalWLeft(elem.W, p);
					double[] derWL = face.evalDerWLeft(elem.W, 0);
					double[] stressVector = eqn.stressVector(WL, derWL, elem.n[faceIdx][p]);
//					double pressure = eqn.pressure(WL);
//					double pressure = eqn.pressure(elem.calculateAvgW());
					for (int d = 0; d < dim; d++) {
//						stressVectors[globalIdx][d] = pressure * elem.n[faceIdx][p][d];
						stressVectors[globalIdx][d] = -stressVector[d];
						force[d] -= stressVector[d] * face.JacobianFace[p] * face.weightsFace[p];
					}
					
					System.arraycopy(face.coordinatesFace[p], 0, stressVectorPositions[globalIdx], 0, dim);
				}
			}
			
			fluidForces[b] = new FluidForces(stressVectors, stressVectorPositions);
			fluidForces[b].force = force;
		}
	}
	
	private BoundaryIndices[] getBoundaryIndexesOfBodies(Element[] elems) {	

//		 number of faces at each body
		int[] nFacesPerBody = new int[nBodies];
		int[] nIntegrPoints = new int[nBodies];
		for (Element elem : elems) {
			for (int r = 0; r < elem.nFaces; r++) {
				if (elem.TEale[r] > 1 && elem.insideMetisDomain) {
					int bodyIdx = elem.TEale[r] - 2;  // indices in TEale start (wierdly) at 2					
					nFacesPerBody[bodyIdx]++;
					nIntegrPoints[bodyIdx] += elem.Int.faces[r].nIntEdge;
				}						
			}
		}
		
		BoundaryIndices[] boundaryInds = new BoundaryIndices[nBodies];
		for (int b = 0; b < nBodies; b++) {
			boundaryInds[b] = new BoundaryIndices(nFacesPerBody[b]);
			boundaryInds[b].nIntegrationPoints = nIntegrPoints[b];
		}
		Arrays.fill(nFacesPerBody, 0);
		
		for (int k = 0; k < elems.length; k++) {
			for (int r = 0; r < elems[k].nFaces; r++) {
				if (elems[k].TEale[r] > 1 && elems[k].insideMetisDomain) {
					int bodyIdx = elems[k].TEale[r] - 2;  // indices in TEale start (wierdly) at 2		
					int i = nFacesPerBody[bodyIdx]++;
					boundaryInds[bodyIdx].elemInds[i] = k;
					boundaryInds[bodyIdx].faceInds[i] = r;
				}						
			}
		}
		
		return boundaryInds;
	}
	
	private class BoundaryIndices {
		int[] elemInds;
		int[] faceInds;
		int nIntegrationPoints;

		public BoundaryIndices(int len) {
			elemInds = new int[len];
			faceInds = new int[len];
		}
	}
	
	public void calculateForces000(Element[] elems, MeshMove[] mshMov) {
		FaceIndexes[][] boundaryIndexesOfBodies = getBoundaryIndexesOfBodies000(elems);
		
		int dim = elems[0].dim;
		for (int b = 0; b < nBodies; b++) {	
			FaceIndexes[] boundaryIndexes = boundaryIndexesOfBodies[b];
			
			double[][] stressVectors = new double[boundaryIndexes.length][dim];
			double[][] stressVectorPositions = new double[boundaryIndexes.length][dim];
			for (int i = 0; i < boundaryIndexes.length; i++) {
				Element elem = elems[boundaryIndexes[i].elemIdx];
				int faceIdx = boundaryIndexes[i].faceIdx;
	//			Element elem = elems[]
	//			stressVector = new
	//			boundaryForce = new double[s][2 * dim]; // [sigmaX, sigmaY, sigmaZ, xes, yes, zes]
	//			faceIndexes = new int[s][3];

//				Integration.Face face = elem.Int.faces[faceIdx];
	//						double[] WL = face.evalWLeft(elem.W, 0);
	//						double[] derWL = face.evalDerWLeft(elem.W, 0);
	//						double[] stressVector = eqn.stressVector(WL, derWL, elem.n[r][0]);
				double pressure = eqn.pressure(elem.calculateAvgW());
				for (int d = 0; d < dim; d++) {
	//							boundaryForce[s][d] = pressure;
					stressVectors[i][d] = pressure * elem.n[faceIdx][0][d];  // prumernou normalu nebo integrovat!!!
	//							boundaryForce[s][d] = -stressVector[d];
				}
				System.arraycopy(elem.Xes[faceIdx], 0, stressVectorPositions[i], 0, dim);
			}
			
			fluidForces[b] = new FluidForces(stressVectors, stressVectorPositions);
		}
	}		
	
	private FaceIndexes[][] getBoundaryIndexesOfBodies000(Element[] elems) {
//		// number of faces at each body
//		int[] nFacesPerBody = new int[nBodies];
//		for (Element elem : elems) {
//			for (int r = 0; r < elem.nFaces; r++) {
//				if (elem.TEale[r] > 1 && elem.insideMetisDomain) {
//					int bodyIdx = elem.TEale[r] - 2;  // indices in TEale start (wierdly) at 2					
//					nFacesPerBody[bodyIdx]++;
//				}						
//			}
//		}
//		
//		// indices of phaces which belong to a specific body (structure)
//		int[][] body2FaceIndexMap = new int[nBodies][];
//		for (int b = 0; b < nBodies; b++) {
//			body2FaceIndexMap[b] = new int[nFacesPerBody[b]];
//		}		

		ArrayList<FaceIndexes>[] faceIndexesAtBodies = new ArrayList[nBodies];
		for (int b = 0; b < nBodies; b++) {
			faceIndexesAtBodies[b] = new ArrayList<>();
		}
		
		for (int k = 0; k < elems.length; k++) {
			for (int r = 0; r < elems[k].nFaces; r++) {
				if (elems[k].TEale[r] > 1 && elems[k].insideMetisDomain) {
					int bodyIdx = elems[k].TEale[r] - 2;  // indices in TEale start (wierdly) at 2					
					faceIndexesAtBodies[bodyIdx].add(new FaceIndexes(k, r));
				}						
			}
		}
		
		FaceIndexes[][] faceIndexesAtBodiesArr = new FaceIndexes[nBodies][];
		for (int b = 0; b < nBodies; b++) {
			faceIndexesAtBodiesArr[b] = faceIndexesAtBodies[b].toArray(new FaceIndexes[0]);
		}
		
		return faceIndexesAtBodiesArr;
	}
	
	private class FaceIndexes {
		int elemIdx;
		int faceIdx;
		int nIntegrPoints;
		
		private FaceIndexes(int elemIdx, int faceIdx) {
			this.elemIdx = elemIdx;
			this.faceIdx = faceIdx;
		}
		
		private FaceIndexes(int elemIdx, int faceIdx, int nIntegrPoints) {
			this(elemIdx, faceIdx);
			this.nIntegrPoints = nIntegrPoints;
		}
	}		
}
