package litempi;

import java.io.*;
import java.net.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MPISlave {

    private static final Logger LOG = LoggerFactory.getLogger(MPISlave.class);

    private Socket socket;
    private ObjectInputStream in;
    private ObjectOutputStream out;

    public MPISlave(String IP, int port) throws MPIException {
        try {
            // if(String.isEmpty(IP)){
            //    sockets[id] = listener.accept();
            // } else {
            //    socket = new Socket();
            //    socket.connect(new InetSocketAddress(IP, port), 1000);
            // {
            //
            socket = new Socket();
            socket.connect(new InetSocketAddress(IP, port), 1000);
            socket.setTcpNoDelay(true);
            socket.setKeepAlive(true);

            /* receive the test message */
            in = new ObjectInputStream(new BufferedInputStream(socket.getInputStream()));
            String testMessage = (String) in.readUnshared();
            LOG.debug("received test message: {}", testMessage);

            /* echo the test message */
            out = new ObjectOutputStream(new BufferedOutputStream(socket.getOutputStream()));
            out.writeUnshared(testMessage);
            out.flush();
            out.reset();

        } catch (IOException | ClassNotFoundException ex) {
            throw new MPIException("error while connecting to " + IP + ":" + port, ex);
        }
    }

    // zasilani zpravy klient
    public synchronized void send(MPIMessage msg) throws MPIException {
        try {
//            ObjectOutputStream out = new ObjectOutputStream(new BufferedOutputStream(socket.getOutputStream()));
            out.writeUnshared(msg);
            out.flush();
        } catch (IOException ex) {
            throw new MPIException("error while sending message", ex);
        }
    }

    // vrati zpravu klient
    public MPIMessage receive() throws MPIException {
        MPIMessage msg = null;
        try {
//            ObjectInputStream in = new ObjectInputStream(new BufferedInputStream(socket.getInputStream()));
//            ObjectInputStream in = new ObjectInputStream(socket.getInputStream());
            msg = (MPIMessage) in.readUnshared();
            //in.reset();
        } catch (IOException | ClassNotFoundException ex) {
            throw new MPIException("error while receiveing message", ex);
        }
        return msg;
    }

    public void reset() throws MPIException {
        try {
            out.flush();
            out.reset();
        } catch (IOException ex) {
            throw new MPIException("error while sending message", ex);
        }
    }

    public void close() throws MPIException {
        try {
            socket.close();
        } catch (Exception ex) {
            throw new MPIException("error while closing socket", ex);
        }
    }
}
