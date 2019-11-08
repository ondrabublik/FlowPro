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
abstract public class Explicit extends TimeIntegrationElement {
    
    Explicit() {
        super();
    }
    
    abstract public int getNumberOfSteps();
    
    abstract public void computeExplicitStep(int step, double dt);
    
    public boolean isImplicit(){
        return false;
    }
    
    public String getLocalSolverType(){
        return "localexplicit";
    }
}
