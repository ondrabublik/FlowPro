package flowpro.core.parallel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
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
    private final List<String> ipList;
    
    public IpAddressContainer(String slaveListFileName, String names2IPMapFileName) throws IOException {         
        loadNames2IPMap(names2IPMapFileName);
        ipList = IpAddressContainer.loadPCList(slaveListFileName, names2IpMap);
    }
    
    public int size() {
        return ipList.size();
    }
    
    public String getIp(int id) {
        return ipList.get(id);
    }
    
    public String getName(int id) {
        String ip = ipList.get(id);
        return ip2NamesMap.getOrDefault(ip, ip);
    }
    
    public String getNameByIp(String ip) {
        return names2IpMap.getOrDefault(ip, ip);
    }
    
    private void loadNames2IPMap(String pcNameFile) throws IOException {
        names2IpMap = new HashMap<>();
        ip2NamesMap = new HashMap<>();
        
        try (BufferedReader reader = new BufferedReader(new FileReader(pcNameFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] token = line.split(" ");
                names2IpMap.put(token[0], token[1]);
                ip2NamesMap.put(token[1], token[0]);
            }
        }
    }

    public static List<String> loadPCList(String ipListFile, Map<String, String> names2IPMap) throws IOException {        
        List<String> ipList = new ArrayList<>();
        
        try (BufferedReader reader = new BufferedReader(new FileReader(ipListFile))) {
            String line;
            InetAddressValidator validator = InetAddressValidator.getInstance();
            while ((line = reader.readLine()) != null) {
                String[] pcNames = line.split(" ");
                for (String pcName : pcNames) {
                    if (names2IPMap.containsKey(pcName)) {
                        ipList.add(names2IPMap.get(pcName));
                    } else {
                        if (validator.isValid(pcName)) {
                            ipList.add(pcName);
                        } else {
                        throw new IOException("error while reading " + ipListFile
                                + ": invalid IP address or unknown PC name \'" + pcName + "\'");
                        }
                    }
                }
            }                        
        }

        return ipList;        
    }
}
