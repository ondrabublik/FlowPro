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
    int nW;
    int nRHS;
    
    TimeIntegration(Parameters par){
        this.method = par.timeMethod;
        this.order = par.orderInTime;
    }
    
    abstract public void init(FlowProProperties props);
    
    abstract public double[] getWcoefficient();
    
    abstract public double[] getRHScoefficient();
    
    abstract public int getWHistoryLength();
    
    abstract public int getRHSHistoryLength();
}
