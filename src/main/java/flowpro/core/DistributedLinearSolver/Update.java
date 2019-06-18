/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.DistributedLinearSolver;

import java.io.Serializable;

class Update implements Serializable {

    double[] y;
    int i;

    Update(double[] y, int i) {
        this.y = y;
        this.i = i;
    }
}
