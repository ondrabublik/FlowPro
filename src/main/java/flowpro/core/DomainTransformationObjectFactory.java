/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core;

import flowpro.api.DomainTransformationObject;
import flowpro.api.FlowProProperties;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLClassLoader;

/**
 *
 * @author obublik
 */
public class DomainTransformationObjectFactory {

    public DomainTransformationObject getDomainTransformationObject(String parameterFilePath, URL[] jarURLList) throws IOException {
        FlowProProperties props = new FlowProProperties();

        try {
            props.load(new FileInputStream(parameterFilePath));
        } catch (IOException ex) {
            throw new IOException("unable to load file " + parameterFilePath, ex);
        }

        String simpleClassName;
        try {
            simpleClassName = props.getString("domainTransformationObject");
        } catch (IOException ex) {
            throw new IOException("file " + parameterFilePath
                    + " has a wrong format: " + ex.getMessage(), ex);
        }

        String className = "flowpro.user.domaintransformationobject." + simpleClassName;
        try {
            Class<DomainTransformationObject> dtoClass = (Class<DomainTransformationObject>) Class.forName(className,true, new URLClassLoader(jarURLList));
            DomainTransformationObject dto = (DomainTransformationObject) dtoClass.newInstance();

            dto.init(props);
            return dto;
        } catch (ClassNotFoundException ex) {
            throw new IOException("class \"" + className + "\" not found, parameter \"domainTransformationObject\" in file "
                    + parameterFilePath + " must have the same name as the class that defines the model", ex);
        } catch (InstantiationException | IllegalAccessException | SecurityException | IllegalArgumentException ex) {
            throw new IOException("error while loading class \"" + className + "\": " + ex, ex);
        }
    }
}
