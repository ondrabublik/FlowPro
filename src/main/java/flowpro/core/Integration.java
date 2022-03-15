package flowpro.core;

import flowpro.api.Mat;
import flowpro.core.transformation.Transformation;
import flowpro.core.transformation.FaceTransformation;
import flowpro.core.quadrature.Quadrature;
import flowpro.core.quadrature.QuadratureCentral;
import flowpro.core.basis.Basis;
import flowpro.core.element.Element;
import flowpro.core.elementType.ElementType;

/**
 *
 * @author obublik
 */
public class Integration {

    int[] TT;
    int[] TP;
    int[] TEshift;
    double[][] shift;
    QuadratureCentral qRules;
    public Quadrature quadVolume;
    public Transformation transform;

    // volume integral
    public int dimension;
    public int nIntVolume;
    public double[] JacobianVolume;
    public double[] weightsVolume;
    public double[][] basisVolume;
    public double[][][] dXibasisVolume;
    public double[][][] dXbasisVolume;
    public double[][] interpolantVolume;

    // edge integral
    int nEdges;
    public Face[] faces;

    public Integration(ElementType elemType, int dimension, Basis basis, Transformation transform, int[] TT, int[] TEshift, double[][] shift, QuadratureCentral qRules) {
        this.transform = transform;
        this.dimension = dimension;
        this.TT = TT;
        this.TEshift = TEshift;
        this.shift = shift;
        this.quadVolume = elemType.getQVolumeRule(qRules);
        this.qRules = qRules;
        nEdges = TT.length;

        // volume
        nIntVolume = quadVolume.nPoints;
        weightsVolume = quadVolume.weights;
        JacobianVolume = transform.jacobian(quadVolume);
        interpolantVolume = transform.getInterpolant(quadVolume);
        basisVolume = basis.getBasisXi(quadVolume.getCoords());
        dXibasisVolume = basis.getDerBasisXi(quadVolume.getCoords(), dimension);
        dXbasisVolume = transform.transformBasis(quadVolume.getCoords(), dXibasisVolume);
        //edge
        faces = new Face[nEdges];
    }

    public void recalculateVolumeIntegrals(Element elem) {
        transform = elem.transform;
        JacobianVolume = transform.jacobian(quadVolume);
        dXbasisVolume = transform.transformBasis(quadVolume.getCoords(), dXibasisVolume);
        for (int i = 0; i < nEdges; i++) {
            faces[i].setFaceTransform(transform);
        }
    }

    public void recalculateFaceIntegrals(Element elem) {
        for (int i = 0; i < nEdges; i++) {
            faces[i].recalculateFaceInt(elem);
        }
    }

    public void initNeighbours(Element[] elems, Element elem) {
        Element elemRight = null;
        for (int k = 0; k < nEdges; k++) {
            if (TT[k] > -1) {
                elemRight = elems[TT[k]];
            }
            // tady by mel vstupovat vetsi z radu sousednich elementu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            faces[k] = new Face(elem.elemType, k, elem.basis, TT[k], TEshift[k], shift, elemRight);
        }
    }

    public class Face {

        public Quadrature quadFace;
        public int nIntEdge;
        private final int TT;
        private int TEshift;
        private double[][] shift;
        private final Element elemRight;
        public int nVerticesEdge;
        public int faceType;
        public int[] faceIndexes;
        public double[] JacobianFace;
        public double[] weightsFace;
        public double[][] coordinatesFace;
        public double[][] coordinatesFaceXiRight;
        public double[][] coordinatesFaceXiRef;
        public double[][] basisFaceLeft, basisFaceRight;
        public double[][][] dXibasisFaceLeft, dXibasisFaceRight;
        
        public double[][][] dXbasisFaceLeft, dXbasisFaceRight;
        public double[][] interpolantFace;
        public FaceTransformation faceTransform;

        Face(ElementType elemType, int k, Basis basis, int TT, int TEshift, double[][] shift, Element elemRight) {
            this.TT = TT;
            this.TEshift = TEshift;
            this.shift = shift;
            this.faceIndexes = elemType.getFaceIndexes(k);
            this.elemRight = elemRight;
            int faceOrder = elemType.faceQuadratureOrder;
            if (elemRight != null) {
                faceOrder = Math.max(elemType.faceQuadratureOrder, elemRight.elemType.faceQuadratureOrder);
            }

            nVerticesEdge = faceIndexes.length;
            faceType = elemType.getFaceType(k);
            quadFace = elemType.getQFaceRule(faceType, qRules, faceOrder);
            faceTransform = elemType.getFaceTransformation(faceType, transform, elemType.getFaceIndexes(k));
            coordinatesFace = faceTransform.getX(quadFace.getCoords());

            nIntEdge = quadFace.nPoints;
            weightsFace = quadFace.weights;
            JacobianFace = faceTransform.jacobian(quadFace);
            coordinatesFaceXiRef = faceTransform.getXiRef(quadFace.getCoords());
            basisFaceLeft = basis.getBasisXi(coordinatesFaceXiRef);
            interpolantFace = faceTransform.getInterpolant(quadFace);
            dXibasisFaceLeft = basis.getDerBasisXi(coordinatesFaceXiRef, dimension);
            dXbasisFaceLeft = transform.transformBasis(coordinatesFaceXiRef, dXibasisFaceLeft);

            if (TT > -1) {
                if (TEshift > 0) { // periodicity shift
                    for (int i = 0; i < coordinatesFace.length; i++) {
                        for (int j = 0; j < coordinatesFace[i].length; j++) {
                            coordinatesFace[i][j] = coordinatesFace[i][j] + shift[TEshift - 1][j];
                        }
                    }
                }
                coordinatesFaceXiRight = elemRight.transform.getXi(coordinatesFace);
                basisFaceRight = elemRight.basis.getBasisXi(coordinatesFaceXiRight);
                dXibasisFaceRight = elemRight.basis.getDerBasisXi(coordinatesFaceXiRight, dimension);
                dXbasisFaceRight = elemRight.transform.transformBasis(coordinatesFaceXiRight, dXibasisFaceRight);
            }
        }

        void recalculateFaceInt(Element elem) {
            coordinatesFace = faceTransform.getX(quadFace.getCoords());
            JacobianFace = faceTransform.jacobian(quadFace);
            dXbasisFaceLeft = transform.transformBasis(coordinatesFaceXiRef, dXibasisFaceLeft);
            if (TT > -1) {
                dXbasisFaceRight = elemRight.transform.transformBasis(coordinatesFaceXiRight, dXibasisFaceRight);
            }
        }

        void setFaceTransform(Transformation transform) {
            faceTransform.setVolumeTransformation(transform);
        }
		
		public double[] evalWLeft(double[] coeffW, int p) {
			double[] base = basisFaceLeft[p];
			int nBasis = base.length;
			int nEqs = coeffW.length / nBasis;
			
			double[] w = new double[nEqs];
			for (int m = 0; m < nEqs; m++) {
				for (int i = 0; i < nBasis; i++) {
					w[m] = w[m] + coeffW[m * nBasis + i] * base[i];
				}
			}		
			return w;
		}

		public double[] evalDerWLeft(double[] coeffW, int p) {
			double[][] derBase = dXbasisFaceLeft[p];
			int nBasis = derBase.length;
			int dim = derBase[0].length;
			int nEqs = coeffW.length / nBasis;
			
			double[] derW = new double[nEqs * dim];
			for (int m = 0; m < nEqs; m++) {
				for (int i = 0; i < nBasis; i++) {
					for (int d = 0; d < dim; ++d) {
						derW[nEqs * d + m] += coeffW[m * nBasis + i] * derBase[i][d];
					}
				}
			}		
			return derW;
		}
		
		public double[] integrateLeft(Integrand integrand, double[] coeffW, double[][] normals) {
			double[] result = new double[dimension];
			for (int p = 0; p < nIntEdge; ++p) {
				double[] w = evalWLeft(coeffW, p);
				double[] dw = evalDerWLeft(coeffW, p);
				
				double[] val = integrand.fun(w, dw, normals[p]);
				for (int d = 0; d < dimension; ++d) {
					result[d] += JacobianFace[p] * weightsFace[p] * val[d];
				}
			}
			return result;
		}
    }
	
	public interface Integrand {
		public double[] fun(double[] w, double[] dw, double[] normal);
	}
}
