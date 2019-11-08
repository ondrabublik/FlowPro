/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import java.util.Arrays;

/**
 *
 * @author obublik
 */
public class RK4 extends Explicit {

    double[] V;
    double[][] RwN;
    double[] W, K1, K2, K3, K4, Wo;
    double[][] iM;

    RK4() {
        super();
    }

    public void init(FlowProProperties props) {
        V = new double[nBasis * nEqs];
        RwN = new double[nFaces][];
        K1 = new double[nBasis * nEqs];
        K2 = new double[nBasis * nEqs];
        K3 = new double[nBasis * nEqs];
        K4 = new double[nBasis * nEqs];
    }

    public int getOrder() {
        return 4;
    }

    public int getNumberOfSteps() {
        return 8;
    }

    public void computeExplicitStep(int step, double dt) {

        W = elem.W;
        Wo = elem.Wo;
        iM = elem.iM;
        elem.limiter(false);
        switch (step) {
            case 1:
                Arrays.fill(K1, 0);
                elem.residuum(V, K1, RwN);
                K1 = Mat.times(iM, K1);
                break;
            case 2:
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt / 2 * K1[j];
                }
                break;
            case 3:
                Arrays.fill(K2, 0);
                elem.residuum(V, K2, RwN);
                K2 = Mat.times(iM, K2);
                break;
            case 4:
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt / 2 * K2[j];
                }
                break;
            case 5:
                Arrays.fill(K3, 0);
                elem.residuum(V, K3, RwN);
                K3 = Mat.times(iM, K3);
                break;
            case 6:
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt * K3[j];
                }
                break;
            case 7:
                Arrays.fill(K4, 0);
                elem.residuum(V, K4, RwN);
                K4 = Mat.times(iM, K4);
                break;
            case 8:
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt / 6 * (K1[j] + 2 * K2[j] + 2 * K3[j] + K4[j]);
                }
                break;
        }
    }
}
