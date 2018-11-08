package flowpro.core.parallel;

//import dgfem2d.Mesh;
//import litempi.MPIMessage;

public class Tag {
    // tags from master to slave
    public static final int CLOSE = 0;
    public static final int INIT_DATA = 1;
    public static final int TIME_STEP_REQ = 2;
    public static final int ASSEMBLE_INTERIOR = 3;
    public static final int ASSEMBLE_OVERLAP_and_SOLVE = 13;
    public static final int ASSEMBLE = 21;
    public static final int SOLVE = 4;
    public static final int LIMITER_and_NEXT_TIME_LEVEL = 50;
    public static final int INTERIOR_LIMITER_and_NEXT_TIME_LEVEL = 5;
    public static final int OVERLAP_UPDATE_and_LIMITER_and_NEXT_TIME = 15;
    public static final int UPDATE_NEWTON = 70;
    public static final int PREVIOUS_TIME_LEVEL = 6;
    public static final int DATA_REQ = 7;
    public static final int DATA_MASTER_TO_SLAVE = 8;
    public static final int SOLUTION_REQ = 9;
    public static final int LIMITERS = 10;
    public static final int ALE_CALCULATE_FORCES = 11;
    public static final int ALE_NEW_MESH_POSITION = 12;
    public static final int ALE_NEXT_TIME_LEVEL = 13;
    public static final int INTERIOR_UPDATE = 14;
    public static final int RESET_OUTPUT_WRITER = 16;
    public static final int GET_MEMORY = 17;
    public static final int GETSOLUTIONMONITOR = 18;
    public static final int SETSOLUTIONMONITOR = 19;
    public static final int GMRES2SLAVE = 20;
    
    // tags from slave to master
    public static final int DATA_INITIALIZED = -1;
    public static final int TIME_STEP = -2;
    public static final int GMRES_RESULT = -4;
    public static final int RESIDUUM = -5;
    public static final int DATA_SLAVE_TO_MASTER = -7;
    public static final int DATA_UPDATED = -8;
    public static final int SOLUTION = -9;
    public static final int FORCES = -10;
    public static final int MEMORY_USAGE = -11;
    public static final int ASSEMBELED = -20;
    public static final int OVERLAP_UPDATED = -21;
    public static final int INTERIOR_UPDATED = -22;
    public static final int MESH_POSITION_UPDATED = -23;
    public static final int NEW_MESH_POSITION_UPDATED = -24;
    public static final int NEWTON_UPDATED = -30;
    public static final int LOCALSOLUTIONMONITOR = -31;
    public static final int SOLUTIONMONITORSET = -32;
    public static final int GMRES2MASTER = -33;
}
