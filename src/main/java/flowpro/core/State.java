package flowpro.core;

import flowpro.api.FlowProProperties;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
public class State implements Serializable {

    public final String stateFilePath;
    public final int order; // order of accuracy in space
    private boolean hasOrderChanged;

    int steps; // number of time steps that has already been taken
    double t;  // physical time
    double cfl;
    double residuum;
    long executionTime;  // current execution time
    long transferTime;   // current transfer time

    long initExecutionTime; // initial execution time
    long initTransferTime;  // initial transfer time

    public State(String stateFilePath, int order) throws IOException {
        this.stateFilePath = stateFilePath;
        this.order = order;
        hasOrderChanged = false;

        steps = 0;
        t = 0.0;
        cfl = 0.0;
        residuum = Double.MAX_VALUE;
        initExecutionTime = 0;
        executionTime = 0;
        initExecutionTime = 0;
        transferTime = 0;
    }

    public void load() throws IOException {
        try {
            FlowProProperties stateProperties = new FlowProProperties();
            stateProperties.load(new FileInputStream(stateFilePath));

            t = stateProperties.getDouble("t");
            steps = stateProperties.getInt("steps");
            initExecutionTime = stateProperties.getLong("CPU");
            initTransferTime = stateProperties.getLong("transfer");
            cfl = stateProperties.getDouble("CFL");

            int oldOrder = stateProperties.getInt("order");
            if (oldOrder != order) {
                hasOrderChanged = true;
            }
        } catch (IOException ex) {
            throw new IOException("file " + stateFilePath + " has a wrong format: "
                    + ex.getMessage());
        }
    }

    public void save() throws IOException {
        FlowProProperties output = new FlowProProperties();

        output.setProperty("t", Double.toString(t));
        output.setProperty("steps", Integer.toString(steps));
        output.setProperty("CPU", Long.toString(getOverallExecutionTime()));
        output.setProperty("transfer", Long.toString(getOverallTransferTime()));
        output.setProperty("CFL", Double.toString(cfl));
        output.setProperty("residuum", Double.toString(residuum));

        output.setProperty("order", Integer.toString(order));

        output.store(new FileOutputStream(stateFilePath), null);
    }

    public boolean hasOrderChanged() {
        return hasOrderChanged;
    }

    public long getOverallExecutionTime() {
        return initExecutionTime + executionTime;
    }

    public long getOverallTransferTime() {
        return initTransferTime + transferTime;
    }
}
