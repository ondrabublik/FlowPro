/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.solver;

/**
 *
 * @author obublik
 */
public class CFLSetup {
    
    double maxCFL;
    boolean varyCFL;
    double res = -1;
    double dres = 0;
    double alfa = 0.1;
    int logResOld;

    CFLSetup(double maxCFL, boolean varyCFL) {
        this.maxCFL = maxCFL;
        this.varyCFL = varyCFL;
    }

    double getCFL(double actualCFL, double residuum) {
        if (varyCFL) {
            if (res == -1) {
                res = residuum;
                logResOld = 1000;
            } else {
                res = alfa * res + (1 - alfa) * residuum; // low pass filter
            }
            int logRes = (int) Math.log(res);
            if (logRes < logResOld) {
                maxCFL *= 1.3;
            } else {
                maxCFL *= 0.8;
            }
            logResOld = logRes;
            actualCFL += maxCFL / 5;
            if (actualCFL > maxCFL) {
                actualCFL = maxCFL;
            }
            return actualCFL;
        } else {
            actualCFL += maxCFL / 20;
            if (actualCFL > maxCFL) {
                actualCFL = maxCFL;
            }
            return actualCFL;
        }
    }

    double reduceCFL(double actualCFL) {
        return actualCFL / 1.5;
    }
}
