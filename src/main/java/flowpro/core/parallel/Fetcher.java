package flowpro.core.parallel;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import javax.net.ssl.SSLSocket;
import javax.net.ssl.SSLSocketFactory;
import org.apache.commons.lang3.tuple.MutablePair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author pecka
 */
public class Fetcher {
    
    private static final Logger LOG = LoggerFactory.getLogger(Fetcher.class);
    
    public static final int TIME_OUT = 1000;

    private final int fetcherPort;
    private final int nSlaves;
    private final Map<String, String> ip2NodeNameMap;
    private final String publicKeyFileName;
    private SSLSocket[] sockets;
    private List<MutablePair<String, Integer>> reachedNodes;
    
    public Fetcher(int nSlaves, Map<String, String> ip2NodeNameMap,
            int fetcherPort, String publicKeyFileName) throws IOException {
        this.nSlaves = nSlaves;
        this.ip2NodeNameMap = ip2NodeNameMap;
        this.fetcherPort = fetcherPort;
        this.publicKeyFileName = publicKeyFileName;        
        
        if (publicKeyFileName == null) {
            throw new IOException("public key was not specified");
        }
    }
    
//    public static void main(String[] args) throws IOException {
////        System.out.println("PC names:");
////        Map<String, String> names2IPMap = Fetcher.loadNames2IPMap("network/pcNames.txt");
////        System.out.println(names2IPMap.toString());
////        System.out.println("PC list:");
////        List<String> pcList = loadPCList("matlab/pclist.txt", names2IPMap);
////        System.out.println(pcList.toString());
//
//        String masterIP = "127.0.0.1";
//        String masterPort = "9001";
//        String parallelSolverType = "test";
//        String arguments = "slave " + masterIP + " " + masterPort + " " + parallelSolverType;              
//        ZipFile zip = new ZipFile("FlowPro.zip", "FlowPro.jar", arguments);
//        IpAddressContainer ipAddresses = new IpAddressContainer("matlab/pclist.txt", "network/pcNames.txt");
//        Fetcher fetcher = new Fetcher(1, ipAddresses, 5555, "testkeystore.ks");
//        fetcher.fetch(zip);
//    }
    
    public List<MutablePair<String, Integer>> establishConnections(List<MutablePair<String, Integer>> nodeList) {
        System.setProperty("javax.net.ssl.trustStore", publicKeyFileName);
        SSLSocketFactory socketFactory = (SSLSocketFactory) SSLSocketFactory.getDefault();
        List<SSLSocket> socketList = new ArrayList<>();
        reachedNodes = new ArrayList<>();
        
        int id = 0;
        int nReachedSlaves = 0;
        while (nReachedSlaves < nSlaves && id < nodeList.size()) {
            String ip = nodeList.get(id).getKey();
            try {                                         
                SSLSocket socket = (SSLSocket) socketFactory.createSocket();
                socket.setTcpNoDelay(true);                        
                socket.connect(new InetSocketAddress(ip, fetcherPort), TIME_OUT);
                socketList.add(socket);
                
                MutablePair<String, Integer> node = nodeList.get(id);
                
                nReachedSlaves += node.getValue();
                
                int nLocalSlaves = Math.min(node.getValue() - (nReachedSlaves - nSlaves), node.getValue());
                reachedNodes.add(new MutablePair<>(node.getKey(), nLocalSlaves));
                
                LOG.info("{}. node {} with {} slave(s) is ready", id+1,
                        ip2NodeNameMap.getOrDefault(ip, ip), reachedNodes.get(id).getValue());
                id++;
            } catch (IOException ex) {
                LOG.info("{} could not be reached: {}", ip2NodeNameMap.getOrDefault(ip, ip), ex.getMessage());                
            }
        }
        
        if (nReachedSlaves < nSlaves) {
            throw new RuntimeException("not enough available slaves");            
        }
        
        sockets = socketList.toArray(new SSLSocket[0]);
        
        return reachedNodes;
    }
    
    private void fetchFile(File zipFile, OutputStream out) throws IOException {
        try (InputStream in = new FileInputStream(zipFile)) {
            byte[] bytes = new byte[16 * 1024];
            int count;
            while ((count = in.read(bytes)) > 0) {
                out.write(bytes, 0, count);
            }
            out.flush();            
        }
    }
    
    public void fetch(File zipFile, AppInfo appInfo) throws IOException {              
        
        try {
            for (int id = 0; id < sockets.length; id++) {
                MutablePair<String, Integer> node = reachedNodes.get(id);
                try (ObjectOutputStream out = new ObjectOutputStream(sockets[id].getOutputStream())) {
                    out.write(1);
                    appInfo.nSlaves = node.getValue();
                    out.writeUnshared(appInfo);
                    fetchFile(zipFile, out);
                    out.flush();
                    LOG.debug("{}/{} zip file has been sent to {}", id+1, nSlaves, node.getKey());
                }
            }
            LOG.info("zip file has been sent to {} slave(s)", nSlaves);
        } catch (IOException ex) {
            throw new IOException("error occured while zip files were being fetched", ex);
        }
    }
}
