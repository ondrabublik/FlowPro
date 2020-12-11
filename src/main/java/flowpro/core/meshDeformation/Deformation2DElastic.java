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
	
	private static final double[][] EYE = new double[][] {{1, 0}, {0, 1}};
	
	private int dim = 2;

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
	
	double[] getStressTensor(Element elem, double[] WL,  double[] derWL) {
		double[][] stressTensor = new double[dim][];
		
		for (int d = 0; d < dim; d++) {
			stressTensor[d] = eqn.stressVector(WL, derWL, EYE[d]);
		}
		
		return new double[] {stressTensor[0][0], stressTensor[1][1], stressTensor[0][1]};
	}
	
	@Override
	public void calculateForces(Element[] elems, MeshMove[] mshMov) {
		BoundaryIndices[] boundaryIndexesOfBodies = getBoundaryIndexesOfBodies(elems);
		
//		int dim = elems[0].dim;
		for (int b = 0; b < nBodies; b++) {	
			BoundaryIndices boundaryIndexes = boundaryIndexesOfBodies[b];
			
			double[] force = new double[dim];
			double[][] stressVectors = new double[boundaryIndexes.nIntegrationPoints][dim];
			double[][] stressTensors = new double[boundaryIndexes.nIntegrationPoints][3];
			double[][] stressVectorPositions = new double[boundaryIndexes.nIntegrationPoints][dim];			
			
			for (int i = 0; i < boundaryIndexes.elemInds.length; i++) {
				Element elem = elems[boundaryIndexes.elemInds[i]];
				int faceIdx = boundaryIndexes.faceInds[i];
				Integration.Face face = elem.Int.faces[faceIdx];
				for (int p = 0; p < face.nIntEdge; p++) {
					int globalIdx = i * face.nIntEdge + p;
					double[] WL = face.evalWLeft(elem.W, p);
					double[] derWL = face.evalDerWLeft(elem.W, p);
					double[] stressVector = eqn.stressVector(WL, derWL, elem.n[faceIdx][p]);
					stressTensors[globalIdx] = getStressTensor(elem, WL, derWL);
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
			fluidForces[b].stressTensors = stressTensors;
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
