package flowpro.core.parallel;

import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.InetSocketAddress;
import javax.net.ssl.SSLSocket;
import javax.net.ssl.SSLSocketFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author pecka
 */
public class Fetcher {
    
    private static final Logger LOG = LoggerFactory.getLogger(Fetcher.class);
    
    public static final int TIME_OUT = 1000;

    private final int fetcherPort = 5555;
    private final int nSlaves;
    
    private final IpAddressContainer ipAddresses;
    
    public Fetcher(int nSlaves, IpAddressContainer ipAddresses, int fetcherPort) throws IOException {
        this.nSlaves = nSlaves;
        this.ipAddresses = ipAddresses;
        
        if (nSlaves > ipAddresses.size()) {
            throw new IOException("IP list contains less PC's than required");
        }
    }
    
    public static void main(String[] args) throws IOException {
//        System.out.println("PC names:");
//        Map<String, String> names2IPMap = Fetcher.loadNames2IPMap("network/pcNames.txt");
//        System.out.println(names2IPMap.toString());
//        System.out.println("PC list:");
//        List<String> pcList = loadPCList("matlab/pclist.txt", names2IPMap);
//        System.out.println(pcList.toString());
        String masterIP = "127.0.0.1";
        String masterPort = "9001";
        String parallelSolverType = "test";
        String arguments = "slave " + masterIP + " " + masterPort + " " + parallelSolverType;              
        ZipFile zip = new ZipFile("FlowPro.zip", "FlowPro.jar", arguments);
        IpAddressContainer ipAddresses = new IpAddressContainer();
        Fetcher fetcher = new Fetcher(1, ipAddresses, 5555);
        fetcher.fetch(zip);
    }
    
    private SSLSocket[] establishConnections() {
        System.setProperty("javax.net.ssl.trustStore", "testkeystore.ks");
        SSLSocketFactory socketFactory = (SSLSocketFactory) SSLSocketFactory.getDefault();
        SSLSocket[] sockets = new SSLSocket[nSlaves];
        
        int id, ipIdx;
        for (id = 0, ipIdx = 0; id < nSlaves && ipIdx < ipAddresses.size(); id++, ipIdx++) {
            try {                                         
            
                SSLSocket socket = (SSLSocket) socketFactory.createSocket();
                socket.setTcpNoDelay(true);                              
                socket.connect(new InetSocketAddress(ipAddresses.getIp(ipIdx), fetcherPort), TIME_OUT);                
                sockets[id] = socket;

                LOG.debug("{}/{} just connected to {}", id+1, nSlaves, ipAddresses.getName(ipIdx));                           
            } catch (IOException ex) {
                LOG.info("{} could not be reached: {}", ipAddresses.getName(ipIdx), ex.getMessage());
                --id;
            }
        }
        
        if (id != nSlaves) {
            throw new RuntimeException("not enough available slaves");            
        }
        
        return sockets;
    }
    
    private void fetchFile(ZipFile zipFile, DataOutputStream out) throws IOException {
        try (InputStream in = new FileInputStream(zipFile.zipFile)) {
            // send name of the .jar file and arguments
            out.writeUTF(zipFile.jarFileName);
            out.writeUTF(zipFile.args);

            // send the actual .zip file
            byte[] bytes = new byte[16 * 1024];
            int count;
            while ((count = in.read(bytes)) > 0) {
                out.write(bytes, 0, count);
            }
            out.flush();            
        }
    }
    
    public void fetch(ZipFile zipFile) throws IOException {
        
        SSLSocket[] sockets = establishConnections();
        LOG.info("connection with {} slave(s) has been successfully established", nSlaves);        
        
        try {
            for (int id = 0; id < nSlaves; id++) {
                try (DataOutputStream out = new DataOutputStream(sockets[id].getOutputStream())) {
                    out.write(1);
                    fetchFile(zipFile, out);
                    LOG.info("{}/{} zip file has been sent to {}", id+1, nSlaves, ipAddresses.getName(id));
                }
            }
            LOG.info("zip file has been sent to {} slave(s)", nSlaves);
        } catch (IOException ex) {
            throw new IOException("error occured while zip files were being fetched", ex);
        }
    }
}
