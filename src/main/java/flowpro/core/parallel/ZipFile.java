package flowpro.core.parallel;

import java.io.File;

public class ZipFile {
    public final File zipFile;
    public final String jarFileName;
    public final String args;

    public ZipFile(File file, String jarFileName, String args) {
        this.zipFile = file;
        this.jarFileName = jarFileName;
        this.args = args;
    }

    public ZipFile(String fileName, String jarFileName, String args) throws IllegalArgumentException { 
        zipFile = new File(fileName);
        this.jarFileName = jarFileName;
        this.args = args;

        if (!zipFile.isFile()) {
            throw new IllegalArgumentException("file " + fileName + " does not exist");
        }
    }

    public String getName() {
        return zipFile.getName();
    }
}
