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
public class RK2Element extends TimeIntegrationElement {
    
    RK2Element(Element elem){
        super(elem);
    }
    
    public void init(){
        
    }
    
    public void computeExplicitStep(double dt) {
        double[] V = new double[nBasis * nEqs];
        double[] Rw = new double[nBasis * nEqs];
        double[][] RwN = new double[nFaces][];
        elem.limiter(false);
        double[] W = elem.W;
        double[] Wo = elem.Wo;
        double[][] iM = elem.iM;
        switch (elem.par.orderInTime) {
            case 2: // RK2
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt / 2 * Rw[j];
                }

                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt * Rw[j];
                }
                break;
            case 3: // RK3
                double[] W1 = new double[W.length];
                double[] W2 = new double[W.length];
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W[j] = Wo[j] + dt * Rw[j];
                    W1[j] = W[j];
                }

                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W[j] = 0.75 * Wo[j] + 0.25 * (W1[j] + dt * Rw[j]);
                }

                Arrays.fill(Rw, 0);
                elem.residuum(V, Rw, RwN);
                Rw = Mat.times(iM, Rw);
                for (int j = 0; j < W.length; j++) {
                    W[j] = 1.0 / 3 * Wo[j] + 2.0 / 3 * (W2[j] + dt * Rw[j]);
                }
        }
        //eqn.limitUnphysicalValues(calculateAvgW(), W, nBasis);
    }
}
