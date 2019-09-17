/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.element;

/**
 *
 * @author obublik
 */
public abstract class TimeIntegration {
    boolean isImplicit;
    Element elem;
    int nEqs;
    int nBasis;
    int nFaces;
    
    TimeIntegration(Element elem){
        this.elem = elem;
        nEqs = elem.nEqs;
        nBasis = elem.nBasis;
        nFaces = elem.nFaces;
    }
    
    abstract public void init();
}