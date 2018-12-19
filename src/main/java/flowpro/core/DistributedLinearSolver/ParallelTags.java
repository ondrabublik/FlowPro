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
public class ParallelTags {
    // tags from master to slave
    public static final int UPDATE_MATRIXES_AND_VECTORS = -1;
    public static final int NORM = 2;
    public static final int SMULT = 0;
    public static final int MULT = 1;
    public static final int FILLV = 8;
    public static final int UPDATE_W = 7;
    public static final int UPDATE_SOLUTION = 6;
    public static final int SAVE_DATA = 3;
    public static final int LOAD_X = 4;
    public static final int LOAD_W = 41;
    public static final int PING = 55;
    
    // tags from slave to master
}
