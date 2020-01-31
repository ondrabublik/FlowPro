package litempi;

import flowpro.core.parallel.IpAddressContainer;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.HashMap;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MPIMaster {

    private static final Logger LOG = LoggerFactory.getLogger(MPIMaster.class);
    
    public static final int TIME_OUT = 30000;

    public final int nSlaves;
    private Socket[] sockets;
    private ObjectOutputStream[] outStreams;
    private ObjectInputStream[] inStreams;
    private final IpAddressContainer ipAddresses;

    // waiting messages from each slaves
    private final Map<Integer, MPIMessage>[] waitingMsgs;

    /**
     * Waits for slaves to connect (old master-server mode).
     * 
     * @param nSlaves
     * @param masterPort
     * @param ipAddressContainer
     * @throws MPIException 
     */
    public MPIMaster(int nSlaves, int masterPort, IpAddressContainer ipAddressContainer) throws MPIException {
        this.nSlaves = nSlaves;
        this.ipAddresses = ipAddressContainer;
        sockets = new Socket[nSlaves];
        outStreams = new ObjectOutputStream[nSlaves];
        inStreams = new ObjectInputStream[nSlaves];

        LOG.info("lite MPI has started and is waiting for slaves to connect...");

        try (ServerSocket listener = new ServerSocket(masterPort)) {
            for (int id = 0; id < nSlaves; id++) {
                sockets[id] = listener.accept();
                sockets[id].setTcpNoDelay(true);
                sockets[id].setKeepAlive(true);

                String testString = sockets[id].getInetAddress().getHostAddress();

                /* sending the test message */
                outStreams[id] = new ObjectOutputStream(new BufferedOutputStream(sockets[id].getOutputStream()));
                outStreams[id].writeUnshared(testString);
                outStreams[id].flush();
                outStreams[id].reset();

                /* receiving echo of the test message */
                inStreams[id] = new ObjectInputStream(new BufferedInputStream(sockets[id].getInputStream()));
                String ipEcho = (String) inStreams[id].readUnshared();

                LOG.info("{}. slave {} is ready", id + 1, ipEcho);
            }
        } catch (IOException | ClassNotFoundException ex) {
            throw new MPIException("error while connecting to slaves", ex);
        }

        waitingMsgs = new Map[nSlaves];
        for (int i = 0; i < nSlaves; ++i) {
            waitingMsgs[i] = new HashMap<>();
        }
    }
    
    /**
     * Connects to slaves (new master-slave mode).
     * 
     * @param ipAddresses
     * @param startPort
     * @throws MPIException 
     */
    public MPIMaster(IpAddressContainer ipAddresses, int startPort) throws MPIException {               
        LOG.info("Lite MPI has started and is about to connect to slaves");
        
        this.ipAddresses = ipAddresses;
        InetSocketAddress[] inetAddresses = ipAddresses.getInetAddresses(startPort);
        this.nSlaves = inetAddresses.length;
        
        try {
            sockets = new Socket[nSlaves];
            outStreams = new ObjectOutputStream[nSlaves];
            inStreams = new ObjectInputStream[nSlaves];
            
            for (int id = 0; id < nSlaves; id++) {
                sockets[id] = new Socket();

                InetSocketAddress addr = inetAddresses[id];
                LOG.debug("master is going to attempt to connect to {} on port {}", addr.getAddress(), addr.getPort());
                sockets[id].connect(new InetSocketAddress(addr.getAddress(), addr.getPort()), TIME_OUT);
                LOG.debug("master has just connected to {}", addr.getAddress());
                sockets[id].setTcpNoDelay(true);
                sockets[id].setKeepAlive(true);

                String testString = sockets[id].getInetAddress().getHostAddress();
                
                LOG.debug("sending test string \'{}\'", testString);
                /* sending the test message */
                outStreams[id] = new ObjectOutputStream(new BufferedOutputStream(sockets[id].getOutputStream()));
                outStreams[id].writeUnshared(testString);
                outStreams[id].flush();
                outStreams[id].reset();

                /* receiving echo of the test message */
                inStreams[id] = new ObjectInputStream(new BufferedInputStream(sockets[id].getInputStream()));
                String testStringEcho = (String) inStreams[id].readUnshared();                
                LOG.debug("received echo of the test string \'{}\'", testStringEcho);
                
                if (!testString.equals(testStringEcho)) {
                    throw new MPIException("test strings do not match (" + testString + " vs. " + testStringEcho + " )");
                }
                
                LOG.info("{}/{} slave {} is ready", id + 1, nSlaves, addr.getAddress());
            }
        } catch (IOException | ClassNotFoundException ex) {
            throw new MPIException("error while connecting to slaves " + ex.getMessage(), ex);
        }

        waitingMsgs = new Map[nSlaves];
        for (int i = 0; i < nSlaves; ++i) {
            waitingMsgs[i] = new HashMap<>();
        }
    }

    public void send(MPIMessage msg, int id) throws MPIException {
        try {
//            LOG.debug("sending {}", msg.tag);
            outStreams[id].writeUnshared(msg);
            outStreams[id].flush();

        } catch (IOException ex) {
            throw new MPIException("error occurred while sending message to "
                    + ipAddresses.getNameByIp(sockets[id].getInetAddress().getHostAddress()), ex);
        }
    }

    class ThreadSend extends Thread {

        int i;
        MPIMessage msg;

        ThreadSend(int i, MPIMessage msg) {
            this.i = i;
            this.msg = msg;
        }

        @Override
        public void run() {
            try {
                outStreams[i].writeUnshared(msg);
                outStreams[i].flush();

            } catch (IOException ex) {
                LOG.error("error occurred while sending message to "
                        + ipAddresses.getNameByIp(sockets[i].getInetAddress().getHostAddress()) + ex);
            }
        }
    }

    public void parallelSendAll(MPIMessage msg) throws MPIException {
        ThreadSend[] ts = new ThreadSend[nSlaves];
        for (int i = 0; i < nSlaves; ++i) {
            ts[i] = new ThreadSend(i, msg);
            ts[i].start();
        }
        try {
            for (int i = 0; i < nSlaves; ++i) {
                ts[i].join();
            }
        } catch (java.lang.InterruptedException ex) {
            LOG.error("parallel");
        }
    }

    public void sendAll(MPIMessage msg) throws MPIException {
        for (int i = 0; i < nSlaves; ++i) {
            try {
//                LOG.debug("sending {}", msg.tag);
                outStreams[i].writeUnshared(msg);
                outStreams[i].flush();

            } catch (IOException ex) {
                throw new MPIException("error occurred while sending message to "
                        + ipAddresses.getNameByIp(sockets[i].getInetAddress().getHostAddress()), ex);
            }
        }
    }

    public void printWaitingMsgs() {
        for (int i = 0; i < nSlaves; ++i) {
            System.out.print(i + "-th slave: ");
            for (Map.Entry<Integer, MPIMessage> entry : waitingMsgs[i].entrySet()) {
                System.out.print(entry.getKey() + " ");
            }
            System.out.println();
        }
    }

    // vrati zpravu server
    public MPIMessage receive(int id, int tag) throws MPIException {
        try {
            MPIMessage msg = waitingMsgs[id].remove(tag);

            if (msg == null) {
                msg = (MPIMessage) inStreams[id].readUnshared();

                while (msg.tag != tag) {
                    waitingMsgs[id].put(msg.tag, msg);
                    msg = (MPIMessage) inStreams[id].readUnshared();
                }
            }

            if (msg instanceof MPIExceptionMessage) {
                MPIExceptionMessage exMsg = (MPIExceptionMessage) msg;
                throw new MPIException("exception thrown by "
                        + ipAddresses.getNameByIp(sockets[id].getInetAddress().getHostAddress()), exMsg.getException());
            }

            return msg;
        } catch (IOException | ClassNotFoundException ex) {
            throw new MPIException("error occurred while receiveing message form "
                    + ipAddresses.getNameByIp(sockets[id].getInetAddress().getHostAddress()), ex);
        }

    }

    public void waitForAll(int tag) throws MPIException {
        for (int i = 0; i < nSlaves; ++i) {
            receive(i, tag);
        }
    }

    public MPIMessage[] receiveAll(int tag) throws MPIException {
        MPIMessage[] msg = new MPIMessage[nSlaves];
        for (int i = 0; i < nSlaves; ++i) {
            msg[i] = receive(i, tag);
        }
        return msg;
    }

    public double[] receiveAllDouble(int tag) throws MPIException {
        double[] array = new double[nSlaves];
        for (int i = 0; i < nSlaves; ++i) {
            array[i] = (double) receive(i, tag).getData();
        }
        return array;
    }

    public double receiveAllDoubleSum(int tag) throws MPIException {
        double sum = 0;
        for (int i = 0; i < nSlaves; ++i) {
            sum += (double) receive(i, tag).getData();
        }
        return sum;
    }

    public double[] receiveAllDoubleArraySum(int tag) throws MPIException {
        double[] sum = null;
        for (int i = 0; i < nSlaves; ++i) {
            if (i == 0) {
                sum = (double[]) receive(i, tag).getData();
            } else {
                double[] array = (double[]) receive(i, tag).getData();
                for (int j = 0; j < sum.length; j++) {
                    sum[j] += array[j];
                }
            }

        }
        return sum;
    }

    public void reset() throws MPIException {
        for (int i = 0; i < nSlaves; ++i) {
            try {
                outStreams[i].flush();
                outStreams[i].reset();
            } catch (IOException ex) {
                throw new MPIException("error occurred while sending message to "
                        + ipAddresses.getNameByIp(sockets[i].getInetAddress().getHostAddress()), ex);
            }
        }
    }

    public void close() throws MPIException {
        try {
            for (Socket soc : sockets) {
                soc.close();
            }
            LOG.info("Lite MPI has shut down");
        } catch (IOException ex) {
            throw new MPIException("error occurred while closing socket", ex);
        }
    }
}
