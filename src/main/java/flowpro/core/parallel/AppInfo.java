package flowpro.core.parallel;

import java.io.Serializable;

public class AppInfo implements Serializable {
    public final String jarFileName;
    public final String[] args;
    public String appName;
    public int nSlaves;

    public AppInfo(String jarFileName, String[] args, String appName) {
        this.jarFileName = jarFileName;
        this.args = args;
        this.appName = appName;
    }
}
