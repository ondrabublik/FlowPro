/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core;

/**
 *
 * @author obublik
 */
public class TimeIntegration {
    String method;
    int order;
    double[][] MWHistory;
    double[][] RHSHistory;
    double[] timeHistory;
    
    TimeIntegration(Parameters par){
        this.method = par.timeMethod;
        this.order = par.orderInTime;
    }
}
