
package flowpro.core.element;

import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class RK2 extends Explicit {

    double[] V, Rw;
    double[][] RwN;
    double[] W1, Wn;
    double[][] iM;

    RK2() {
        super();
    }

    public void init(FlowProProperties props) {
        V = new double[nBasis * nEqs];
        Rw = new double[nBasis * nEqs];
        RwN = new double[nFaces][];
        W1 = new double[nBasis * nEqs];
        Wn = new double[nBasis * nEqs];
    }

    public int getOrder(){
        return 2;
    }
    
    public int getNumberOfSteps(){
        return 4;
    }
    
    public void computeExplicitStep(int step, double dt) {

        double[] W = elem.W;
        double[] Wo = elem.Wo;
        iM = elem.iM;
        switch (step) {
            case 1:
                elem.limiter(false);
                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W1[j] = Wo[j] + dt/2 * Rw[j];
                }
                break;
            case 2:
                System.arraycopy(W1, 0, W, 0, W.length);
                break;
            case 3:
                elem.limiter(false);
                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    Wn[j] = Wo[j] + dt * Rw[j];
                }
                break;
            case 4:
                System.arraycopy(Wn, 0, W, 0, W.length);
                break;
        }
        //eqn.limitUnphysicalValues(calculateAvgW(), W, nBasis);
    }
}