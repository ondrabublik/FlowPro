/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.DistributedLinearSolver;

/**
 *
 * @author obublik
 */
public class Watch {
    private long startTime;
    
    Watch(){
        
    }
    
    void start(){
        startTime = System.nanoTime();
    }
    
    void stop(){
        System.out.println("Time: " + (System.nanoTime() - startTime)/1000 + "us");
    }
}
