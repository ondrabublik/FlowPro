package flowpro.core.meshDeformation;

import flowpro.api.Equation;
import flowpro.api.FluidForces;
import flowpro.api.Mat;
import flowpro.api.MeshMove;
import flowpro.core.Integration;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.HashSet;
import java.util.Iterator;

/**
 *
 * @author obublik
 */
public class Deformation2DElastic1 extends Deformation {

	public double[][] boundaryForce;  // smazat
	public int[][] faceIndexes;  // smazat
	
	private final FluidForces[] fluidForces;
	
	private static final double[][] EYE = new double[][] {{1, 0}, {0, 1}};
	
	private static final int dim = 2;
	
	private boolean isInitialized;
	private double[][] staticBndPoints;
	private double[][] bndPoints;
	private double[][] bndDisplacement;
	
	double[][] displacement;
//	double[][] displacementOld;
	
	double[][] points;
	
	private int step;
//	private static double p0 = 285.7142857142857;  // CFD2 and FSI2
//	private static double p0 = 1.143e+06;  // CFD3 and FSI3
//	private static double p0 = 0;

	public Deformation2DElastic1(Parameters par, Equation eqn, double[][] points, int[][] TP, int[][] TEale) throws IOException {
		super(par, eqn, TEale);
		
		this.points = points;
		
		isInitialized = false;
		step = 0;
//		bndPoints = new double[nBodies][][];
//		bndDisplacement = new double[nBodies][][];
//		displacement = new double[nBodies][][];
//		displacementOld = new double[nBodies][][];
		
//		Set<Integer>[] vertexIndicesOnBodies = new Set[nBodies];
//		for (int b = 0; b < nBodies; b++) {
//			vertexIndicesOnBodies[b] = new HashSet<>();
//		}		
//		for (int k = 0; k < TEale.length; k++) {  // element index
//			for (int i = 0; i < TEale[k].length; i++) {  // face index
//				if (TEale[k][i] > 0) {
//					int bodyIdx = TEale[k][i] - 1;
//					vertexIndicesOnBodies[bodyIdx].add(TP[k][i]);
//				}
//			}
//		}
		
		Set<Integer> staticBndNodes = new HashSet<>();
		Set<Integer> movingBndNodes = new HashSet<>();
		for (int k = 0; k < TEale.length; k++) {  // element index
			for (int i = 0; i < TEale[k].length; i++) {  // face index
				if (TEale[k][i] == 1 || TEale[k][i] == 3) {  // TEale[k][i] == 3 jen docasne pro Turek-Hron benchmark (staticky valec)
					staticBndNodes.add(TP[k][i]);
					int iNext = (i + 1) % TP[k].length;
					staticBndNodes.add(TP[k][iNext]);
				} else if (TEale[k][i] > 1) {
					movingBndNodes.add(TP[k][i]);
					int iNext = (i + 1) % TP[k].length;
					movingBndNodes.add(TP[k][iNext]);
				}
//				if (TEale[k][i] > 0) {
//					int bodyIdx = TEale[k][i] - 1;
//					vertexIndicesOnBodies[bodyIdx].add(TP[k][i]);
//				}
			}
		}
		staticBndNodes.removeAll(movingBndNodes);
		
		staticBndPoints = new double[staticBndNodes.size()][dim];
		Iterator<Integer> iterator = staticBndNodes.iterator();
		for (int i = 0; iterator.hasNext(); i++) {
			int elemIdx = iterator.next();
			for (int d = 0; d < dim; d++) {
				staticBndPoints[i][d] = points[elemIdx][d];
			}
		}		
		
		fluidForces = new FluidForces[nBodies];
	}

	private double[][] computeRadialBasisFunctionCoefficients(double[][] points, double[][] b) {
		int nPoints = points.length;
		double[][] rbfCoeff = new double[dim][nPoints];
		double[][] A = new double[nPoints][nPoints];		
		for (int i = 0; i < nPoints; i++) {
			for (int j = 0; j < nPoints; j++) {
				A[i][j] = radialBasisFunction(points[i], points[j]);
			}
		}
		
		int n = b[0].length;
		double[] b0 = new double[n];
		for (int d = 0; d < dim; d++) {
			double[][] A0 = Mat.copyMat(A);			
            System.arraycopy(b[d], 0, b0, 0, n);
//			rbfCoeff[d] = Mat.cg(A, b[d], 1e-6, 500);
			rbfCoeff[d] = Mat.lsolve(A0, b0, b0.length);
		}
		
		return rbfCoeff;
	}
	
	private double[][] evaluateRadialBasisFunction(double[][] rbfCoeffs, double[][] oldPoints, double[][] newPoints) {		
		double[][] newValues =  new double[newPoints.length][dim];
		
		for (int p = 0; p < rbfCoeffs[0].length; p++) {
			for (int j = 0; j < newPoints.length; j++) {
				double rbf = radialBasisFunction(newPoints[j], oldPoints[p]);
				for (int d = 0; d < dim; d++) {
					newValues[j][d] += rbfCoeffs[d][p] * rbf;
				}
			}
		}
		
		return newValues;
	}
	
//@Override
	public void newMeshPositionAndVelocity1(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {
		
		if (!isInitialized) {
			isInitialized = true;
			
//			for (int b = 0; b < nBodies; b++) {
			int b = 0;
			
			double[][] movingBndPoints = mshMov[b].getBoundaryPoints();

			int len = movingBndPoints.length + staticBndPoints.length;
			bndPoints = new double[len][dim];
			bndDisplacement = new double[dim][len];

			for (int p = 0; p < staticBndPoints.length; p++) {
				System.arraycopy(staticBndPoints[p], 0, bndPoints[p+movingBndPoints.length], 0, dim);
			}

			double[][] dis = mshMov[b].getBoundaryDisplacement();
			displacement = new double[dis.length][dis[0].length];
//			}			
		}

		double a1 = 1.0 / dt;
		double a2 = -1.0 / dt;
		double a3 = 0;

		if (timeOrder > 1) {
			a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
			a2 = -(dt + dto) / (dt * dto);  // -2/dt;
			a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
		}
		
//		for (int b = 0; b < nBodies; b++) {
		int b = 0;
		
		double[][] displacementOld = displacement;
		displacement = mshMov[b].getBoundaryDisplacement();
		double[][] diffBoundaryDisplacement = Mat.minus(displacement, displacementOld);
		double[][] movingBndPoints = mshMov[b].getBoundaryPoints();

		for (int p = 0; p < movingBndPoints.length; p++) {
			for (int d = 0; d < dim; d++) {
				bndPoints[p][d] = movingBndPoints[p][d] + displacement[d][p];
			}
		}

		for (int d = 0; d < dim; d++) {
			System.arraycopy(diffBoundaryDisplacement[d], 0, bndDisplacement[d], 0, diffBoundaryDisplacement[d].length);
		}

		double[][] rbfCoeffs = computeRadialBasisFunctionCoefficients(bndPoints, bndDisplacement);

		for (Element elem : elems) {
			double[][] diffDispacement = evaluateRadialBasisFunction(rbfCoeffs, bndPoints, elem.vertices);
			for (int j = 0; j < elem.nVertices; j++) {
				for (int d = 0; d < dim; d++) {
					elem.vertices[j][d] += diffDispacement[j][d];
				}
			}
		}
//		}

		for (Element elem : elems) {
			for (int j = 0; j < elem.nVertices; j++) {
				for (int d = 0; d < elem.dim; d++) {
					elem.U[j][d] = a1 * elem.vertices[j][d] + a2 * elem.verticesOld[j][d] + a3 * elem.verticesOld2[j][d];
				}
			}
		}
		
		step++;
	}
	
	@Override
	public void newMeshPositionAndVelocity(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {
		
		if (!isInitialized) {
			isInitialized = true;
			
			int movingBndLen = 0;
//			for (int b = 0; b < nBodies; b++) {
			int b = 0;
				double[][] movingBndPoints = mshMov[b].getBoundaryPoints();
				movingBndLen += movingBndPoints.length;
//			}
			int totalLen = movingBndLen + staticBndPoints.length;
			bndPoints = new double[totalLen][dim];
			bndDisplacement = new double[dim][totalLen];

			for (int p = 0; p < staticBndPoints.length; p++) {
				System.arraycopy(staticBndPoints[p], 0, bndPoints[p+movingBndLen], 0, dim);
			}			
		}

		double a1 = 1.0 / dt;
		double a2 = -1.0 / dt;
		double a3 = 0;

		if (timeOrder > 1) {
			a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
			a2 = -(dt + dto) / (dt * dto);  // -2/dt;
			a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
		}
		
		// multiple blending function
//		for (Element elem : elems) {
//			for (int j = 0; j < elem.nVertices; j++) {
//				for (int d = 0; d < elem.dim; d++) {
//					elem.vertices[j][d] = 0;
//				}
//			}
//		}
		
		int pos = 0;
//		for (int b = 0; b < nBodies; b++) {
		int b = 0;
			displacement = mshMov[b].getBoundaryDisplacement();
			double[][] diffBoundaryDisplacement = displacement;
			double[][] movingBndPoints = mshMov[b].getBoundaryPoints();

			for (int p = 0; p < movingBndPoints.length; p++) {
				System.arraycopy(movingBndPoints[p], 0, bndPoints[p], 0, dim);
			}

			for (int d = 0; d < dim; d++) {
				System.arraycopy(diffBoundaryDisplacement[d], 0, bndDisplacement[d], pos, diffBoundaryDisplacement[d].length);
			}
			
			pos += diffBoundaryDisplacement[0].length;
//		}

		double[][] rbfCoeffs = computeRadialBasisFunctionCoefficients(bndPoints, bndDisplacement);

		for (Element elem : elems) {
			double[][] diffDispacement = evaluateRadialBasisFunction(rbfCoeffs, bndPoints, elem.vertices0);
			for (int j = 0; j < elem.nVertices; j++) {
				for (int d = 0; d < dim; d++) {
					elem.vertices[j][d] = elem.vertices0[j][d] + diffDispacement[j][d];
				}
			}
		}
		
//		}

		for (Element elem : elems) {
			for (int j = 0; j < elem.nVertices; j++) {
				for (int d = 0; d < elem.dim; d++) {
					elem.U[j][d] = a1 * elem.vertices[j][d] + a2 * elem.verticesOld[j][d] + a3 * elem.verticesOld2[j][d];
				}
			}
		}
		
		step++;
	}

//	@Override
//	public void newMeshPositionAndVelocity2(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {
//		
//		if (!isInitialized) {
//			for (int b = 0; b < nBodies; b++) {
//				double[][] points = mshMov[b].getBoundaryPoints();
//				
//				int len = staticBndPoints.length + points.length;
//				bndPoints[b] = new double[len][dim];
//				bndDisplacement[b] = new double[dim][len];
//				
//				for (int p = 0; p < points.length; p++) {
//					System.arraycopy(points[p], 0, bndPoints[b][p], 0, dim);
//				}
//				
//				for (int p = 0; p < staticBndPoints.length; p++) {
//					System.arraycopy(staticBndPoints[p], 0, bndPoints[b][p+points.length], 0, dim);
//				}
//			}
//			
//			isInitialized = true;
//		}
//
//		double a1 = 1.0 / dt;
//		double a2 = -1.0 / dt;
//		double a3 = 0;
//
//		if (timeOrder > 1) {
//			a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
//			a2 = -(dt + dto) / (dt * dto);  // -2/dt;
//			a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
//		}
//		
//		// multiple blending function
//		for (Element elem : elems) {
//			for (int j = 0; j < elem.nVertices; j++) {
//				for (int d = 0; d < elem.dim; d++) {
//					elem.vertices[j][d] = 0;
//				}
//			}
//		}
//		
//		for (int b = 0; b < nBodies; b++) {
//			double[][] displacement = mshMov[b].getBoundaryDisplacement();
//			
//			for (int d = 0; d < dim; d++) {
//				System.arraycopy(displacement[d], 0, bndDisplacement[b][d], 0, displacement[d].length);
//			}
//			
//			double[][] rbfCoeffs = computeRadialBasisFunctionCoefficients(bndPoints[b], bndDisplacement[b]);
//			
//			for (Element elem : elems) {		
//				double[][] dispacement = evaluateRadialBasisFunction(rbfCoeffs, bndPoints[b], elem.vertices0);
//				for (int j = 0; j < elem.nVertices; j++) {
//					for (int d = 0; d < dim; d++) {
//						elem.vertices[j][d] += elem.vertices0[j][d] + dispacement[j][d];
//					}
//				}
//			}
//		}
//
//		for (Element elem : elems) {
//			for (int j = 0; j < elem.nVertices; j++) {
//				for (int d = 0; d < elem.dim; d++) {
//					elem.U[j][d] = a1 * elem.vertices[j][d] + a2 * elem.verticesOld[j][d] + a3 * elem.verticesOld2[j][d];
//				}
//			}
//		}
//	}
	
//	@Override
	public void newMeshPositionAndVelocityOld(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {

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
			for (int j = 0; j < elem.nVertices; j++) {
				for (int d = 0; d < elem.dim; d++) {
					elem.vertices[j][d] = 0;
				}
			}
		}
		
		for (int b = 0; b < nBodies; b++) {
			double[][] bndPoints = mshMov[b].getBoundaryPoints();
			double[][] bndDisplacement = mshMov[b].getBoundaryDisplacement();
			double[][] rbfCoeffs = computeRadialBasisFunctionCoefficients(bndPoints, bndDisplacement);
			
			for (Element elem : elems) {		
				double[][] dispacement = evaluateRadialBasisFunction(rbfCoeffs, bndPoints, elem.vertices0);
				for (int j = 0; j < elem.nVertices; j++) {
					for (int d = 0; d < dim; d++) {
						double xNew = elem.vertices0[j][d] + dispacement[j][d] * nBodies;
						elem.vertices[j][d] += ((1 - elem.blendFun[j][b]) * elem.vertices0[j][d] + elem.blendFun[j][b] * xNew) / nBodies;
					}
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
		return 1 - norm * norm * norm;
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
//					double[] stressVector = eqn.stressVector(WL, derWL, elem.n[faceIdx][p]);					
					
					stressTensors[globalIdx] = getStressTensor(elem, WL, derWL);
					
					double[] stressVector = new double[dim];  // stressTensor * normal
					stressVector[0] = stressTensors[globalIdx][0] * elem.n[faceIdx][p][0]
							+ stressTensors[globalIdx][2] * elem.n[faceIdx][p][1];
					stressVector[1] = stressTensors[globalIdx][2] * elem.n[faceIdx][p][0]
							+ stressTensors[globalIdx][1] * elem.n[faceIdx][p][1];
					
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
		
		// jen docasne - je treba pri nacitani rozlisit tuhy, elasticky i stacionarni telesa
		if (nBodies > 1) {
			fluidForces[1].bodyType = FluidForces.BodyType.STATIC;
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
