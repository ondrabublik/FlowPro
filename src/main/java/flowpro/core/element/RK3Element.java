/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.Mat;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class RK3Element extends TimeIntegrationElement {

    double[] V, Rw;
    double[][] RwN;
    double[] W, W1, W2, Wo, Wn;
    double[][] iM;

    RK3Element(Element elem) {
        super(elem);

        V = new double[nBasis * nEqs];
        Rw = new double[nBasis * nEqs];
        RwN = new double[nFaces][];
        W1 = new double[nBasis * nEqs];
        W2 = new double[nBasis * nEqs];
        Wn = new double[nBasis * nEqs];
    }

    public void init() {

    }

    public void computeExplicitStep(int step, double dt) {

        W = elem.W;
        Wo = elem.Wo;
        iM = elem.iM;
        elem.limiter(false);
        switch (step) {
            case 1: // RK3
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W1[j] = Wo[j] + dt * Rw[j];
                }
                break;
            case 11: // RK3
                System.arraycopy(W1, 0, W, 0, W.length);
                break;
            case 2:
                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W2[j] = 0.75 * Wo[j] + 0.25 * (W1[j] + dt * Rw[j]);
                }
                break;
            case 21: // RK3
                System.arraycopy(W2, 0, W, 0, W.length);
                break;
            case 3:
                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    Wn[j] = 1.0 / 3 * Wo[j] + 2.0 / 3 * (W2[j] + dt * Rw[j]);
                }
                break;
            case 31: // RK3
                System.arraycopy(Wn, 0, W, 0, W.length);
                break;
        }
        //eqn.limitUnphysicalValues(calculateAvgW(), W, nBasis);
    }
}
