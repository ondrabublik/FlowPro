package flowpro.core.parallel;

import java.io.Serializable;

public class AppInfo implements Serializable {
    public final String jarFileName;
    public final String[] args;
    public int nSlaves;

    public AppInfo(String jarFileName, String[] args) {
        this.jarFileName = jarFileName;
        this.args = args;
        nSlaves = args.length;
    }
}
