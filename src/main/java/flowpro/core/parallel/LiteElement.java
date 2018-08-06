package flowpro.core.parallel;

/**
 *
 * @author obublik
 */
public class LiteElement implements java.io.Serializable {

    public int index;
    public double[] y;

    public LiteElement() {
    }
    
    public LiteElement(int index, double[] y) {
        this.index = index;
        this.y = y;
    }
}
