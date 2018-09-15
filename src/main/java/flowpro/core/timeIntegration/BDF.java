/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.timeIntegration;

import flowpro.api.FlowProProperties;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class BDF extends TimeIntegration {

    public void init(FlowProProperties props) throws IOException {
        if (props.containsKey("orderInTime")) {
            order = props.getInt("orderInTime");
        } else {
            order = 1;
            System.out.println(" order in time was set to 1! ");
        }

        nMWcoef = order + 1;
        nRHScoef = 0;
        Wcoef = new double[nMWcoef];
        RHScoef = null;
    }

    public void computeTimecoefficient(double[] dt) {
        switch (order) {
            case 1:
                Wcoef[0] = 1 / dt[0];
                Wcoef[1] = -1 / dt[0];
                break;
            case 2:
                Wcoef[0] = (2 * dt[0] + dt[1]) / (dt[0] * (dt[0] + dt[1]));  // 3/(2*dt);
                Wcoef[1] = -(dt[0] + dt[1]) / (dt[0] * dt[1]);  // -2/dt;
                Wcoef[2] = dt[0] / (dt[1] * (dt[0] + dt[1]));  // 1/(2*dt);
                break;
            default:
                throw new RuntimeException("solver supports only first and second order in time");
        }
    }
}
