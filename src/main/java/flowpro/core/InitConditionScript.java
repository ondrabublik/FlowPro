package flowpro.core;

import java.io.FileNotFoundException;
import javax.script.Invocable;
import javax.script.ScriptEngineManager;
import javax.script.ScriptEngine;
import javax.script.ScriptException;

/**
 *
 * @author obublik
 */
public class InitConditionScript {

    private double[][] initW = null;

    InitConditionScript(double[][] PXY, int[][] TP, int nEqs, String scriptLocation) throws FileNotFoundException, ScriptException, NoSuchMethodException {
        
        double[][] init = new double[TP.length][nEqs];
        
        ScriptEngineManager factory = new ScriptEngineManager();
        ScriptEngine engine = factory.getEngineByName("js");
        engine.eval("load(\"" + scriptLocation + "\");");
        Invocable invocable = (Invocable) engine;
        
        int dim = PXY[0].length;
        for (int i = 0; i < TP.length; i++) {
            double[] x = new double[dim];
            int m = TP[i].length;
            for (int j = 0; j < m; j++) {
                for (int d = 0; d < dim; d++) {
                    x[d] += PXY[TP[i][j]][d] / m;
                }
            }
            String[] ls = invocable.invokeFunction("init", x).toString().split(" ");
            for (int k = 0; k < nEqs; k++) {
                init[i][k] = Double.parseDouble(ls[k]);
            }
        }
        initW = init;
    }

    public double[][] getInitW() {
        return initW;
    }
}
