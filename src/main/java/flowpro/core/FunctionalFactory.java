package flowpro.core;

import flowpro.api.FlowProProperties;
import flowpro.api.Functional;
import java.io.FileInputStream;
import java.io.IOException;

/**
 *
 * @author obublik
 */
public class FunctionalFactory {

    public Functional getFunctional(String parameterFilePath) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("functionalOfInterest");
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        }

        String className = "flowpro.user.functional." + simpleClassName;

        try {
            Class<Functional> funClass = (Class<Functional>) Class.forName(className);
//            Constructor constructor = eqnClass.getConstructor();  // new Class[]{FlowProProperties.class}
//            Object obj = constructor.newInstance(props);
            Functional fun = (Functional) funClass.newInstance();

            try {
                fun.init(props);
                return fun;
            } catch (IOException ex) {
                throw new IOException("file " + parameterFilePath
                        + " has a wrong format: " + ex.getMessage());
            }
        } catch (ClassNotFoundException ex) {
            throw new IOException("class \"" + className + "\" not found, parameter \"model\" in file "
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
