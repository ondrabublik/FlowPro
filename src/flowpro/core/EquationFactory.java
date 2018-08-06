package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.api.Equation;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLClassLoader;
//import java.lang.reflect.Constructor;
//import java.lang.reflect.InvocationTargetException;

/**
 *
 * @author ales
 */
public class EquationFactory {

    public Equation getEquation(String parameterFilePath, URL[] jarURLList) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("model");
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        }

        String className = "flowpro.user.equation." + simpleClassName;

        try {
            Class<Equation> eqnClass = (Class<Equation>) Class.forName(className,true, new URLClassLoader(jarURLList));
            Equation eqn = (Equation) eqnClass.newInstance();

            try {
                eqn.init(props);
                return eqn;
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
