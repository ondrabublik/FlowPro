
package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.api.SolutionMonitor;
import java.io.FileInputStream;
import java.io.IOException;
//import java.lang.reflect.Constructor;
//import java.lang.reflect.InvocationTargetException;

/**
 *
 * @author ales
 */
public class SolutionMonitorFactory {

    public SolutionMonitor getSolutionMonitor(String parameterFilePath) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("solutionMonitor").split("\'")[1];
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        } catch (ArrayIndexOutOfBoundsException ex) {
            throw new IOException("file " + parameterFilePath + " has a wrong format: "
                    + "the name of the model is required to be surrounded by apostrophe, e.g. \'myModel\'", ex);
        }

        String className = "flowpro.user.solutionMonitor." + simpleClassName;

        try {
            Class<SolutionMonitor> eqnClass = (Class<SolutionMonitor>) Class.forName(className);
//            Constructor constructor = eqnClass.getConstructor();  // new Class[]{FlowProProperties.class}
//            Object obj = constructor.newInstance(props);
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
//        } catch (NoSuchMethodException ex) {
//            throw new IOException("incorrect implementation of class \"" + className +
//                    "\", one-parameter constructor of form \"public "
//            + simpleClassName + "(" + FlowProProperties.class.getSimpleName() + ")\" is required");
        } catch (InstantiationException | IllegalAccessException | SecurityException | IllegalArgumentException ex) {
//                | InvocationTargetException | NoSuchMethodException ex) {
            throw new IOException("error while loading class \"" + className + "\": " + ex, ex);
        }
    }
}
