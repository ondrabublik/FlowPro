package litempi;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.HashMap;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MPIMaster {

    private static final Logger LOG = LoggerFactory.getLogger(MPIMaster.class);

    public final int nSlaves;
    private Socket[] sockets;
    private ObjectOutputStream[] outStreams;
    private ObjectInputStream[] inStreams;
    
    // waiting messages from each slaves
    private final Map<Integer, MPIMessage>[] waitingMsgs;

    // navazani spojeni s pocitaci
    public MPIMaster(int nSlaves, String masterIP, int masterPort) throws MPIException {
        this.nSlaves = nSlaves;
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

    public void send(MPIMessage msg, int id) throws MPIException {
        try {
//            LOG.debug("sending {}", msg.tag);
            outStreams[id].writeUnshared(msg);
            outStreams[id].flush();

        } catch (IOException ex) {
            throw new MPIException("error while sending message to "
                    + sockets[id].getInetAddress().getHostAddress(), ex);
        }
    }

    public void sendAll(MPIMessage msg) throws MPIException {
        for (int i = 0; i < nSlaves; ++i) {
            try {
//                LOG.debug("sending {}", msg.tag);
                outStreams[i].writeUnshared(msg);
                outStreams[i].flush();

            } catch (IOException ex) {
                throw new MPIException("error while sending message to "
                        + sockets[i].getInetAddress().getHostAddress(), ex);
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
//                LOG.debug("received {}", msg.tag);

                while (msg.tag != tag) {
                    waitingMsgs[id].put(msg.tag, msg);
//                    LOG.debug("added " + msg.tag);
                    msg = (MPIMessage) inStreams[id].readUnshared();
//                    LOG.debug("received {}", msg.tag);
                }
            }

            if (msg instanceof MPIExceptionMessage) {
                MPIExceptionMessage exMsg = (MPIExceptionMessage) msg;
                throw new MPIException("exception thrown by "
                        + sockets[id].getInetAddress().getHostAddress(), exMsg.getException());
            }
            
            return msg;
        } catch (IOException | ClassNotFoundException ex) {
            throw new MPIException("error while receiveing message form "
                    + sockets[id].getInetAddress().getHostAddress(), ex);
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

    public void reset() throws MPIException {
        for (int i = 0; i < nSlaves; ++i) {
            try {
                outStreams[i].flush();
                outStreams[i].reset();
            } catch (IOException ex) {
                throw new MPIException("error while sending message to "
                        + sockets[i].getInetAddress().getHostAddress(), ex);
            }
        }
    }

    public void close() throws MPIException {
        try {
            for (Socket soc : sockets) {
                soc.close();
            }
            LOG.info("lite MPI has shut down");
        } catch (IOException ex) {
            throw new MPIException("error while closing socket", ex);
        }
    }
}
