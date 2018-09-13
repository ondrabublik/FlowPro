/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.timeIntegration;

import flowpro.api.FlowProProperties;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLClassLoader;

/**
 *
 * @author ales
 */
public class TimeIntegrationFactory {

    public TimeIntegration getTimeIntegrationMethod(String parameterFilePath, URL[] jarURLList) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("timeIntegrationMethod");
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        }

        String className = "flowpro.core.timeIntegration" + simpleClassName;

        try {
            Class<TimeIntegration> timClass = (Class<TimeIntegration>) Class.forName(className,true, new URLClassLoader(jarURLList));
            TimeIntegration tim = (TimeIntegration) timClass.newInstance();

//            try {
                tim.init(props);
                return tim;
//            } catch (IOException ex) {
//                throw new IOException("file " + parameterFilePath
//                        + " has a wrong format: " + ex.getMessage());
//            }
        } catch (ClassNotFoundException ex) {
            throw new IOException("class \"" + className + "\" not found, parameter \"model\" in file "
                    + parameterFilePath + " must have the same name as the class that defines the model", ex);
        } catch (InstantiationException | IllegalAccessException | SecurityException | IllegalArgumentException ex) {
            throw new IOException("error while loading class \"" + className + "\": " + ex, ex);
        }
    }
}
