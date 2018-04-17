package flowpro.core.LinearSolvers;

import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;
import java.io.*;
import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author obublik
 */
public class ExternSolver extends LinearSolver {

    int dofs;
    int nnz;
    public double[] IA, JA, HA, b;
    String path = "toolbox/matlabInterface/linearSolver/";

    //Matlab
    MatlabClient mc;

    ExternSolver(Element[] elems, int dofs, Parameters par) {
        this.elems = elems;
        this.dofs = dofs;
        nnz = getNNZ();
        IA = new double[nnz];
        JA = new double[nnz];
        HA = new double[nnz];
        b = new double[dofs];

        // Launching Matlab
        try {
            mc = new MatlabClient();
            mc.init();
        } catch (Exception e) {
            System.out.println("Matlab init error " + e);
        }
    }

    @Override
    public boolean solve(double[] x) {
        try {
            buildMatrix();
            mc.solve(IA, JA, HA, b, x);
            return true;
        } catch (Exception ex) {
            Logger.getLogger(ExternSolver.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }
    }

    public int getNNZ() {
        int s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            s += n * n;
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    s += n * ne;
                }
            }
        }
        return s;
    }

    public void buildMatrix() throws IOException {
        int s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            int[] glob = elem.gi_U;
            double[][] Ad = elem.ADiag;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    IA[s] = glob[i];
                    JA[s] = glob[j];
                    HA[s] = Ad[j][i];
                    s++;
                }
            }
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TT[k] > -1) {
                    int ne = elem.getNEqs() * elems[elem.TT[k]].nBasis;
                    int[] globe = elems[elem.TT[k]].gi_U;
                    double[][] An = elem.ANeighs[k].MR;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < ne; j++) {
                            IA[s] = glob[i];
                            JA[s] = globe[j];
                            HA[s] = An[j][i];
                            s++;
                        }
                    }
                }
            }
        }
        
        s = 0;
        for (Element elem : elems) {
            int n = elem.getNEqs() * elem.nBasis;
            for (int i = 0; i < n; i++) {
                b[s] = elem.RHS_loc[i];
                s++;
            }
        }
    }
}

class MatlabClient {

    Socket socket = null;
    ObjectOutputStream out;
    ObjectInputStream in;

    MatlabClient() {
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

    void solve(double[] IA, double[] JA, double[] HA, double[] b, double[] x) throws IOException, ClassNotFoundException {
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
