/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import java.io.Serializable;

class GMRESdata implements Serializable {
    // initialize workspace

    int m, n;
    double[][] V;
    double[] w, r, aux;

    GMRESdata(int m, int n) {
        w = new double[n];
        r = new double[n];
        aux = new double[n];
        V = new double[m + 1][n];
    }
}
