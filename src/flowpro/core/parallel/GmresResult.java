package flowpro.core.parallel;

import java.io.Serializable;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author ales
 */
public class GmresResult implements Serializable {

    public boolean converges;
    public double schwarzResid;

    public GmresResult(boolean converges, double schwarzResid) {
        this.converges = converges;
        this.schwarzResid = schwarzResid;
    }
}
