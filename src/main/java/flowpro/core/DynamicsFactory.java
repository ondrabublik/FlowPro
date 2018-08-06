package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.api.Dynamics;
import flowpro.api.Equation;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLClassLoader;

/**
 *
 * @author obublik
 */
public class DynamicsFactory {
   public Dynamics getDynamicsModel(String parameterFilePath, URL[] jarURLList, int nBodies, String simulationPath, String geometryPath, Equation eqn) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("dynamicsModel");
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        }

        String className = "flowpro.user.dynamics." + simpleClassName;

        try {
            Class<Dynamics> eqnClass = (Class<Dynamics>) Class.forName(className,true, new URLClassLoader(jarURLList));
            Dynamics dyn = (Dynamics) eqnClass.newInstance();

            try {
                dyn.init(nBodies, simulationPath, geometryPath, eqn);
                return dyn;
            } catch (IOException ex) {
                throw new IOException("file " + parameterFilePath
                        + " has a wrong format: " + ex.getMessage());
            }
        } catch (ClassNotFoundException ex) {
            throw new IOException("class \"" + className + "\" not found, parameter \"model\" in file "
                    + parameterFilePath + " must have the same name as the class that defines the model", ex);
        } catch (InstantiationException | IllegalAccessException | SecurityException | IllegalArgumentException ex) {
            throw new IOException("error while loading class \"" + className + "\": " + ex, ex);
        }
    } 
}