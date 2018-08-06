
package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.api.SolutionMonitor;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLClassLoader;

/**
 *
 * @author ales
 */
public class SolutionMonitorFactory {

    public SolutionMonitor getSolutionMonitor(String parameterFilePath, URL[] jarURLList) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("solutionMonitor");
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        }

        String className = "flowpro.user.solutionMonitor." + simpleClassName;

        try {
            Class<SolutionMonitor> eqnClass = (Class<SolutionMonitor>) Class.forName(className,true, new URLClassLoader(jarURLList));
            SolutionMonitor sMon = (SolutionMonitor) eqnClass.newInstance();

            try {
                sMon.init(props);
                return sMon;
            } catch (IOException ex) {
                throw new IOException("file " + parameterFilePath
                        + " has a wrong format: " + ex.getMessage());
            }
        } catch (ClassNotFoundException ex) {
            throw new IOException("class \"" + className + "\" not found, parameter \"solutionMonitor\" in file "
                    + parameterFilePath + " must have the same name as the class that defines the model", ex);
        } catch (InstantiationException | IllegalAccessException | SecurityException | IllegalArgumentException ex) {
            throw new IOException("error while loading class \"" + className + "\": " + ex, ex);
        }
    }
}
