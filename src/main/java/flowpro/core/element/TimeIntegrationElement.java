/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

import flowpro.api.FlowProProperties;

/**
 *
 * @author obublik
 */
public abstract class TimeIntegrationElement {

    boolean isImplicit;
    Element elem;
    int nEqs;
    int nBasis;
    int nFaces;

    TimeIntegrationElement() {
    }

    void set(Element elem) {
        this.elem = elem;
        nEqs = elem.nEqs;
        nBasis = elem.nBasis;
        nFaces = elem.nFaces;
    }

    abstract public void init(FlowProProperties props);

    abstract public int getOrder();

    abstract public boolean isImplicit();
    
    abstract public String getLocalSolverType();
}
