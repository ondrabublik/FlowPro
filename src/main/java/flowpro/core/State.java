package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.core.element.Element;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
public class State implements Serializable {
	
	private boolean hasOrderChanged;
	
    public final String stateFilePath;
    public final int order; // order of accuracy in space    

	public final boolean isTimeStepFixed;
	public final double dt0;
	private double maxCFL;
	private double dt;
	
    public int steps; // number of time steps that has already been taken
    public double t;  // physical time		
    public double currentCFL;	
    public double residuum;
    public long executionTime;  // current execution time
    public long transferTime;   // current transfer time

    public long initExecutionTime; // initial execution time
    public long initTransferTime;  // initial transfer time
	
    private boolean varyCFL;
    private double res = -1;
    private static final double ALPHA = 0.1;
    private int logResOld;

    public State(String stateFilePath, Parameters par) throws IOException {
        this.stateFilePath = stateFilePath;
        this.order = par.order;
        hasOrderChanged = false;

        steps = 0;
        t = 0.0;
        currentCFL = 0.0;
        residuum = Double.MAX_VALUE;
        initExecutionTime = 0;
        executionTime = 0;
        transferTime = 0;
		
		maxCFL = par.cfl;
		isTimeStepFixed = par.isTimeStepFixed;
		dt0 = par.dt;
		
		if (isTimeStepFixed) {
			dt = dt0;
		}
    }

    public void load() throws IOException {
        try {
            FlowProProperties stateProperties = new FlowProProperties();
            stateProperties.load(new FileInputStream(stateFilePath));

            t = stateProperties.getDouble("t");
            steps = stateProperties.getInt("steps");
            initExecutionTime = stateProperties.getLong("CPU");
            initTransferTime = stateProperties.getLong("transfer");
			
            currentCFL = stateProperties.getDouble("CFL");

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
        output.setProperty("residuum", Double.toString(residuum));
		output.setProperty("CFL", Double.toString(currentCFL));
		output.setProperty("dt", Double.toString(dt));

        output.setProperty("order", Integer.toString(order));

        output.store(new FileOutputStream(stateFilePath), null);
    }
	
	public double getTimeStep(Element[] elems) {
		double temp = Double.MAX_VALUE;
		for (Element elem : elems) {
			double localEigenvalue = elem.calculateMaxEigenvalue();
			temp = Double.min(temp, elem.elemSize / localEigenvalue);
		}
		
		if (isTimeStepFixed) {
			currentCFL = dt / temp;
		} else {
			dt = currentCFL * temp;
		}

        return dt;
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

    public void updateCFL() {
		if (isTimeStepFixed) {
//			if (dt < dt0) {
//				dt *= Double.min(2 * dt, dt0);
//				steps--;
//			}
		} else {
			if (varyCFL) {
				if (res == -1) {
					res = residuum;
					logResOld = 1000;
				} else {
					res = ALPHA * res + (1 - ALPHA) * residuum; // low pass filter
				}
				int logRes = (int) Math.log(res);
				if (logRes < logResOld) {
					maxCFL *= 1.3;
				} else {
					maxCFL *= 0.8;
				}
				logResOld = logRes;
				currentCFL += maxCFL / 5;
				if (currentCFL > maxCFL) {
					currentCFL = maxCFL;
				}
			} else {
				currentCFL += maxCFL / 20;
				if (currentCFL > maxCFL) {
					currentCFL = maxCFL;
				}
			}
		}
    }

    public void reduceCFL() {
		if (isTimeStepFixed) {
//			dt /= 2;
		} else {
			currentCFL /= 1.5;
		}
    }
}
