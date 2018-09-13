/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.timeIntegration;

import flowpro.api.FlowProProperties;
import flowpro.core.Parameters;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
abstract public class TimeIntegration implements Serializable {
    String method;
    int order;
    double[][] WHistory;
    double[][] MWHistory;
    double[][] RHSHistory;
    double[] timeHistory;
    
    TimeIntegration(Parameters par){
        this.method = par.timeMethod;
        this.order = par.orderInTime;
    }
    
    abstract void init(FlowProProperties props);
    
    abstract void saveTimeIteration(double[] W, double[] MW, double[] RHS);
    
    abstract double[] getWnm1();
    
    abstract void getBack();
}
