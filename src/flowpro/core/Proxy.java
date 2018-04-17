package flowpro.core;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.ConnectException;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.net.SocketAddress;

import javax.net.ssl.SSLSocketFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author ales
 */
public class Proxy implements Runnable {

    private static final Logger log = LoggerFactory.getLogger(Proxy.class);
    private final String jarFileName;
    private final String outputDirName;
    private final Object lock;
    private boolean stop;
    private final String ip;
    private final int port;

    public Proxy(String ip, int port, String jarFileName, String outputDirName, Object lock) {
        this.ip = ip;
        this.port = port;
        this.jarFileName = jarFileName;
        this.outputDirName = outputDirName;
        this.lock = lock;
        this.stop = false;
    }

    public void stop() {
        stop = true;
    }

    private static Socket connect(String address, int port) throws IOException {
        SocketAddress addr = new InetSocketAddress(address, port);

        System.setProperty("javax.net.ssl.trustStore", "testkeystore.ks");
        SSLSocketFactory socketFactory = (SSLSocketFactory) SSLSocketFactory.getDefault();
        Socket socket;
        long connectDelay = 1000 + (int)(Math.random() * 3000); ;
        for (int i = 0; i < 20; ++i) {
            try {
                socket = socketFactory.createSocket();
                socket.setTcpNoDelay(true);
                socket.connect(addr, 50);
                return socket;

            } catch (ConnectException ex) {
                log.warn("failed to connect to " + address + ":" + port
                        + ", going to try again in " + connectDelay / 1000. + " seconds");
                try {
                    Thread.sleep(connectDelay);
                } catch (InterruptedException ex2) {
                    System.err.println(ex2);
                }
            }
        }

        throw new IOException("failed to connect to " + address + ":" + port);
    }

//    public static void sendExecutionTime(long executionTime, String ip, int port) {
//        try (Socket socket = connect(ip, port);
//                DataOutputStream out = new DataOutputStream(socket.getOutputStream())) {
//            
//            out.writeLong(executionTime);
//
//        } catch (IOException ex) {
//            System.err.println("Proxy: error while sending data: " + ex);
//        }
//    }

    private void fetchFile(File file, DataOutputStream out) throws IOException {
        try (InputStream in = new FileInputStream(file)) {
            // send name of the file and arguments
            out.writeUTF(file.getName());
            out.writeLong(file.length());

            // send the actual file
            byte[] bytes = new byte[16 * 1024];
            int count;
            while ((count = in.read(bytes)) > 0) {
                out.write(bytes, 0, count);
            }
            out.flush();
        }
    }

    private void fetchDirectoryContent(String dirName, DataOutputStream out) throws IOException {
        File dir = new File(dirName);

        for (File file : dir.listFiles()) {
            fetchFile(file, out);
        }
    }

    @Override
    public void run() {
        try (Socket socket = connect(ip, port);
                DataOutputStream out = new DataOutputStream(socket.getOutputStream())) {

            log.info("Proxy:  just connected to " + socket.getInetAddress().getHostAddress());

            out.writeUTF(jarFileName + " " + outputDirName);

            synchronized (lock) {
                try {
                    lock.wait();
                } catch (InterruptedException ex) {
                    log.error("", ex);
                }
                while (!stop) {
                    try {
                        log.info("Proxy: sending data to " + ip);
                        fetchDirectoryContent(outputDirName, out);
                        lock.wait();
                    } catch (IOException ex) {
                        log.error("Proxy: error while sending data to " + ip + " " + ex.getMessage());
                        log.info("killing computation");
                        System.exit(1);
                    } catch (InterruptedException ex) {
                        log.error("", ex);
                        log.info("killing computation");
                        System.exit(1);
                    }
                }
            }
            //deleteDir(new File(outputDirName));
        } catch (IOException ioexception) {
            log.error((new StringBuilder()).append("Proxy: ").append(ioexception.getMessage()).toString());
        }
    }

    public static void deleteDir(File dir) {
        File[] files = dir.listFiles();
        if (files != null) {  // some JVMs return null for empty dirs
            for (File f : files) {
                if (f.isDirectory()) {
                    deleteDir(f);
                } else {
                    f.delete();
                }
            }
        }
        dir.delete();
    }
}
