/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.timeIntegration;

import flowpro.api.FlowProProperties;
import java.io.IOException;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
abstract public class TimeIntegration implements Serializable {

    int order;
    int nMWcoef;
    int nRHScoef;
    double[] Wcoef;
    double[] RHScoef;

    abstract public void init(FlowProProperties props) throws IOException;

    abstract public void computeTimecoefficient(double[] dt);

    public double[] getWcoefficient() {
        return Wcoef;
    }

    public double[] getRHScoefficient() {
        return RHScoef;
    }

    public int getWHistoryLength() {
        return nMWcoef - 1;
    }

    public int getRHSHistoryLength() {
        return nRHScoef - 1;
    }
}
