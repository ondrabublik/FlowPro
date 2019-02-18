package flowpro.core.elementType;

import flowpro.core.Parameters;
import flowpro.core.basis.*;
import flowpro.core.transformation.*;
import flowpro.core.curvedBoundary.FaceCurvature;
import flowpro.core.quadrature.Quadrature;
import flowpro.core.quadrature.QuadratureCentral;
import java.io.IOException;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
public abstract class ElementType implements Serializable {

    public int nVertices;
    public int nFaces;
    public int order;

    public abstract int numberOfPoints();

    public abstract int numberOfEdges();

    public abstract Basis getBasis(Transformation transform) throws IOException;

    public abstract Quadrature getQVolumeRule(QuadratureCentral qRules);

    public abstract Transformation getVolumeTransformation(double[][] vertices, FaceCurvature fCurv, Parameters par);

    public abstract int getFaceType(int k);

    public abstract int[] getFaceIndexes(int k);

    public FaceTransformation getFaceTransformation(int elemFaceType, Transformation transform, int[] faceIndexes) {
        switch (elemFaceType) {
            case 1:
                return new faceTransformation1Dpoint(transform, faceIndexes);
            case 2:
                return new faceTransformation2Dline(transform, faceIndexes);
            case 3:
                return new faceTransformation3Dtriangle(transform, faceIndexes);
            case 4:
                return new faceTransformation3Dsquare(transform, faceIndexes);
            default:
                throw new RuntimeException("unknown face type");
        }
    }
    
    public Quadrature getQFaceRule(int elemFaceType, QuadratureCentral qRules, int faceOrder) {
        switch (elemFaceType) {
            case 1:
                return qRules.qPoint[faceOrder];
            case 2:
                return qRules.qLine[faceOrder];
            case 3:
                return qRules.qTriangle[faceOrder];
            case 4:
                return qRules.qSquare[faceOrder];
            default:
                throw new RuntimeException("unknown face type");
        }
    }
    
    public static ElementType elementTypeFactory(int type, int order) {
        ElementType elemType = null;
        switch (type) {
            case 1:
                elemType = new pointElement(order);
                break;
            case 2:
                elemType = new lineElement(order);
                break;
            case 3:
                elemType = new triangleElement(order);
                break;
            case 33:
                elemType = new triangleElementUserDef(order);
                break;
            case 31:
                elemType = new curvedBoundaryTriangleElement(order);
                break;
            case 4:
                elemType = new squareElement(order);
                break;
            case 42:
                elemType = new squareElementTaylor(order);
                break;    
            case 5:
                elemType = new tetrahedralElement(order);
                break;
            case 6:
                elemType = new hexahedralElement(order);
                break;
            case 63:
                elemType = new hexahedralElementUserDef(order);
                break; 
            case 7:
                elemType = new prismElement(order);
                break;
        }
        return elemType;
    }

    public static int firstDigit(int number) {
        int digit = number;
        while (number > 0) {
            digit = number;
            number = number / 10;
        }
        return digit;
    }

    public static int[] firstDigit(int[] number) {
        for (int i = 0; i < number.length; i++) {
            number[i] = firstDigit(number[i]);
        }
        return number;
    }
}
