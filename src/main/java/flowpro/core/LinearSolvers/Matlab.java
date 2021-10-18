/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.Parameters;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author obublik
 */
public class Matlab extends LinearSolver {

    int dofs;
    double[] b;

    SparseMatrix A;
    
    //Matlab
    MatlabClient3 mc;

    Matlab(SparseMatrix A, Parameters par) throws IOException {
        //this.dofs = dofs;
        this.A = A;
        b = new double[dofs];

        // Launching Matlab
        try {
            mc = new MatlabClient3();
            mc.init();
        } catch (Exception e) {
            System.out.println("Matlab init error " + e);
        }
    }

    @Override
    public boolean solve(double[] x, double[] b) {
        try {
            A.updateData();
            A.updateB(b);
            mc.solve(A.Icoo, A.Jcoo, A.getData(), b, x);
            return true;
        } catch (Exception ex) {
            Logger.getLogger(Matlab.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }
    }
}

class MatlabClient3 {

    Socket socket = null;
    ObjectOutputStream out;
    ObjectInputStream in;

    MatlabClient3() {
        try (ServerSocket listener = new ServerSocket(5767)) {
            socket = listener.accept();
            socket.setTcpNoDelay(true);
            socket.setKeepAlive(true);
            out = new ObjectOutputStream(new BufferedOutputStream(socket.getOutputStream()));
            out.writeObject("test");
            out.flush();
            in = new ObjectInputStream(new BufferedInputStream(socket.getInputStream()));
            in.readObject();
            System.out.println("Succesfully connect with Matlab ...");
        } catch (Exception e) {
            System.out.println(e);
        }
    }

    void init() throws IOException, ClassNotFoundException {
        out.writeObject("init");
        out.flush();
        in.readObject();
    }

    void solve(int[] IA, int[] JA, double[] HA, double[] b, double[] x) throws IOException, ClassNotFoundException {
        out.writeObject("solve");
        out.flush();
        out.writeUnshared(IA);
        out.writeUnshared(JA);
        out.writeUnshared(HA);
        out.writeUnshared(b);
        out.flush();
        out.reset();
        double[] x2 = (double[]) in.readObject();
        System.arraycopy(x2, 0, x, 0, x.length);
        in.readObject();
    }

    void close() throws IOException {
        socket.close();
    }
}

