package flowpro.core.parallel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.validator.routines.InetAddressValidator;

/**
 *
 * @author pecka
 */
public class IpAddressContainer {

    private Map<String, String> names2IpMap;
    private Map<String, String> ip2NamesMap;
    private final List<Node> nodeList;

    public IpAddressContainer(String slaveListFileName, String names2IPMapFileName) throws IOException {
        loadNames2IPMap(names2IPMapFileName);
        nodeList = IpAddressContainer.loadPCList(slaveListFileName, names2IpMap);
    }
    
    public int maxSlavesPerNode() {
        Comparator<Node> comparator = (n1, n2) -> {
            if (n1.nSlaves == n2.nSlaves) return 0;
            else if (n1.nSlaves == n2.nSlaves) return 1;
            else return -1;
        };
        return Collections.max(nodeList, comparator).nSlaves;
    }

    public int size() {
        return nodeList.size();
    }

    public String getIp(int id) {
        return nodeList.get(id).ip;
    }
    
    public int getNumberOfSlaves(int id) {
        return nodeList.get(id).nSlaves;
    }

    public void removeIp(int id) {
        nodeList.remove(id);
    }

    public String getName(int id) {
        String ip = nodeList.get(id).ip;
        return ip2NamesMap.getOrDefault(ip, ip);
    }

    public String getNameByIp(String ip) {
        return ip2NamesMap.getOrDefault(ip, ip);
    }
    
    public InetSocketAddress[] getInetAddresses(int startPort) {
        List<InetSocketAddress> inetList = new ArrayList<>();
        nodeList.forEach((n) -> {
            for (int i = 0; i < n.nSlaves; i++) {
                inetList.add(new InetSocketAddress(n.ip, startPort+i));
            }
        });
        
        return inetList.toArray(new InetSocketAddress[0]);
    }

    private void loadNames2IPMap(String pcNameFile) throws IOException {
        names2IpMap = new HashMap<>();
        ip2NamesMap = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(pcNameFile))) {
            String line;
            for (int lineIdx = 0; (line = reader.readLine()) != null; lineIdx++) {
                String[] token = line.split(" ");
                
                switch (token.length) {
                    case 0:
                        continue;
                    case 2:
                        names2IpMap.put(token[0], token[1]);
                        ip2NamesMap.put(token[1], token[0]);
                        break;
                    default:
                        throw new IOException("error in file " + pcNameFile + " on line " + lineIdx + 1);
                }                                                  
            }
        }
    }

    private static List<Node> loadPCList(String ipListFile, Map<String, String> names2IPMap) throws IOException {
        List<Node> ipList = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(ipListFile))) {
            String line;
            InetAddressValidator validator = InetAddressValidator.getInstance();
            for (int lineIdx = 0; (line = reader.readLine()) != null; lineIdx++) {
                String[] words = line.trim().split("\\s+");

                int nSlavesPerNode;
                String pcName;
                switch (words.length) {
                    case 0:
                        continue;
                    case 1:
                        nSlavesPerNode = 1;
                        pcName = words[0];
                        break;
                    case 2:                        
                        nSlavesPerNode = Integer.parseInt(words[0]);
                        pcName = words[1];                        
                        break;
                    default:
                        throw new IOException("error in file " + ipListFile + " on line " + (lineIdx + 1));
                }
                
                String ip;
                if (names2IPMap.containsKey(pcName)) {
                    ip = names2IPMap.get(pcName);

                } else if (validator.isValid(pcName)) {
                    ip = pcName;
                } else {
                    throw new IOException("error while reading " + ipListFile
                            + ": invalid IP address or unknown PC name \'" + pcName + "\'");
                }
                ipList.add(new Node(ip, nSlavesPerNode));
            }
        }

        return ipList;
    }
}

class Node {

    final String ip;
    final int nSlaves;

    Node(String ip, int nSlaves) {
        this.ip = ip;
        this.nSlaves = nSlaves;
    }
}
