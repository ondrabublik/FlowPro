package flowpro.core.parallel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.validator.routines.InetAddressValidator;

/**
 *
 * @author pecka
 */
public class IpAddressReader {

    private Map<String, String> nodeName2IpMap;
    private Map<String, String> ip2NodeNameMap;
    private final List<MutablePair<String, Integer>> nodeList;

    public IpAddressReader(String slaveListFileName, String names2IPMapFileName) throws IOException {
        loadNames2IPMap(names2IPMapFileName);
        nodeList = loadPCList(slaveListFileName, nodeName2IpMap);
    }
    
    public List<MutablePair<String, Integer>> getNodeList() {
        return nodeList;
    }
    
    public Map<String, String> getIp2NodeNameMap() {
        return ip2NodeNameMap;
    }    

    private void loadNames2IPMap(String pcNameFile) throws IOException {
        nodeName2IpMap = new HashMap<>();
        ip2NodeNameMap = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(pcNameFile))) {
            String line;
            for (int lineIdx = 0; (line = reader.readLine()) != null; lineIdx++) {
                String[] token = line.split(" ");
                
                switch (token.length) {
                    case 0:
                        continue;
                    case 2:
                        nodeName2IpMap.put(token[0], token[1]);
                        ip2NodeNameMap.put(token[1], token[0]);
                        break;
                    default:
                        throw new IOException("error in file " + pcNameFile + " on line " + lineIdx + 1);
                }                                                  
            }
        }
    }

    private List<MutablePair<String, Integer>> loadPCList(String ipListFile, Map<String, String> names2IPMap) throws IOException {
        List<MutablePair<String, Integer>> ipList = new ArrayList<>();

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
                ipList.add(new MutablePair(ip, nSlavesPerNode));
            }
        }

        return ipList;
    }
}
