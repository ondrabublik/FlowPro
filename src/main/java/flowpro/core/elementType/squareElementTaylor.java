package flowpro.core.elementType;

import flowpro.core.basis.*;
import flowpro.core.transformation.*;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class squareElementTaylor extends squareElement{
    
    public squareElementTaylor(int order, int volumeQuardatureOrder, int faceQuardatureOrder) {
        super(order, volumeQuardatureOrder, faceQuardatureOrder);
    }

    public Basis getBasis(Transformation transform) throws IOException {
        Basis basis = new basis2DsquareTaylor(order);
        basis.calculateCoefficients();
        return basis;
    }
}
