/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.DistributedLinearSolver;

import flowpro.core.Parameters;
import flowpro.core.parallel.Domain;
import flowpro.core.parallel.LiteElement;
import flowpro.core.parallel.Tag;
import java.util.logging.Level;
import java.util.logging.Logger;
import litempi.*;

/**
 *
 * @author obublik
 */
public class ParallelGmresMaster {

    MPIMaster mpi;
    Domain domain;
    LiteElement[] liteElems;

    int n, m, iterationLimit, nThreads;
    double tol;
    double[][] H;
    double[] cs, sn, e1;

    public ParallelGmresMaster(Parameters par, int m, MPIMaster mpi, Domain domain) {
        this.mpi = mpi;
        this.domain = domain;
        iterationLimit = 5;
        tol = par.iterativeSolverTol;
        nThreads = par.nThreads;
        liteElems = new LiteElement[domain.nElems];
        this.m = m;

        // initialize workspace
        H = new double[m + 1][m];
        cs = new double[m];
        sn = new double[m];
        e1 = new double[m + 1];
    }

    public boolean solve() throws MPIException {
        double temp;
        double norm_w = 0;
        double[] s, y;
        double[] Hk = null;
        int gsIter = 1;

        double norm_b = updatePreconditioner();   //----
        if (norm_b == 0) {
            norm_b = 1;
        }
        
        double norm_r = M_bmAx(); //----
        double error = norm_r / norm_b;
        if (error < tol) {
            return true;
        }
        e1[0] = 1;
        for (int iter = 0; iter < iterationLimit; iter++) {          // begin iteration
            if (iter > 0) {
                norm_r = norm('r'); //----
            }
            fillV(0, norm_r); //----
            s = timesVectorByScalar(e1, norm_r);
            for (int i = 0; i < m; i++) {                        // construct orthonormal
                for (int gs = 0; gs < gsIter; gs++) { // gram-schmidt ortogonalisation                 
                    Hk = M_AV(i); //----
                    norm_w = updateW(Hk); //----                 
                }
                for (int k = 0; k <= i; k++) {
                    H[k][i] = Hk[k];
                }
                H[i + 1][i] = norm_w;
                fillV(i + 1, H[i + 1][i]); //----
                for (int k = 0; k <= i - 1; k++) {                              // apply Givens rotation
                    temp = cs[k] * H[k][i] + sn[k] * H[k + 1][i];
                    H[k + 1][i] = -sn[k] * H[k][i] + cs[k] * H[k + 1][i];
                    H[k][i] = temp;
                }
                double[] rot = rotmat(H[i][i], H[i + 1][i]); // form i-th rotation matrix
                cs[i] = rot[0];
                sn[i] = rot[1];
                temp = cs[i] * s[i];                            // approximate residual norm
                s[i + 1] = -sn[i] * s[i];
                s[i] = temp;
                H[i][i] = cs[i] * H[i][i] + sn[i] * H[i + 1][i];
                H[i + 1][i] = 0;
                error = Math.abs(s[i + 1]) / norm_b;
                if (error <= tol) {                        // Update approximation
                    y = lsolve(H, s, i + 1);                 // and exit 
                    updateSolution(y, i); //----
                    break;
                }
            }
            
            if (error <= tol) {
                break;
            }
            y = lsolve(H, s, m);
            updateSolution(y, m - 1); //----
            norm_r = M_bmAx(); //----
            s[m] = norm_r;
            error = s[m] / norm_b;                     // check convergence
            if (error <= tol) {
                break;
            }
        }
        return error <= tol;
    }

    double M_bmAx() throws MPIException {
        dataExchange(-1);
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.SMULT));
        double normR = mpi.receiveAllDoubleSum(Tag.GMRES2MASTER);
        return Math.sqrt(normR);
    }

    double[] M_AV(int i) throws MPIException {
        dataExchange(i);
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.MULT, i));
        return mpi.receiveAllDoubleArraySum(Tag.GMRES2MASTER);
    }

    void fillV(int i, double norm) throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.FILLV, new double[]{(double) i, norm}));
        mpi.waitForAll(Tag.GMRES2MASTER);
    }

    double updateW(double[] H) throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.UPDATE_W, H));
        double normW = mpi.receiveAllDoubleSum(Tag.GMRES2MASTER);
        return Math.sqrt(normW);
    }

    double updatePreconditioner() throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.UPDATE_MATRIXES_AND_VECTORS));
        double normB = mpi.receiveAllDoubleSum(Tag.GMRES2MASTER);
        return Math.sqrt(normB);
    }

    void updateSolution(double[] y, int m) throws MPIException {
        Update up = new Update(y, m);
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.UPDATE_SOLUTION, up));
        mpi.waitForAll(Tag.GMRES2MASTER);
    }

    void dataExchange(int index) throws MPIException {
        boolean parallel = true;

        // downloading data from the slave nodes into the central structure
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.SAVE_DATA, index));
        if (!parallel) {// sequential
            for (int d = 0; d < domain.nDoms; ++d) {
                int[] mapL2G = domain.getSubdomain(d).mapL2G;
                LiteElement[] dataRcv = (LiteElement[]) mpi.receive(d, Tag.GMRES2MASTER).getData();
                for (LiteElement dataRcv1 : dataRcv) {
                    //liteElems[mapL2G[dataRcv1.index]] = new LiteElement(dataRcv1.index, dataRcv1.y);
                    liteElems[mapL2G[dataRcv1.index]] = dataRcv1;
                }
            }
        } else {  // parallel
            ParallelSaving[] ps = new ParallelSaving[domain.nDoms];
            for (int d = 0; d < domain.nDoms; d++) {
                ps[d] = new ParallelSaving(d);
                ps[d].start();
            }
            try {
                for (int d = 0; d < domain.nDoms; d++) {
                    ps[d].join();
                }
            } catch (java.lang.InterruptedException e) {
                System.out.println(e);
            }
        }

        // uploading data from the central structure onto the slave nodes
        if (!parallel) {// sequential
            for (int d = 0; d < domain.nDoms; ++d) {
                Domain.Subdomain subdom = domain.getSubdomain(d);
                int[] mapG2L = subdom.mapG2L;
                int[] load = subdom.load1;
                LiteElement[] dataSend = new LiteElement[load.length];
                for (int i = 0; i < load.length; i++) {
                    //dataSend[i] = new LiteElement(mapG2L[load[i]], liteElems[load[i]].y);
                    dataSend[i] = liteElems[load[i]];
                    dataSend[i].index = mapG2L[load[i]];
                }
                if (index == -1) {
                    mpi.send(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.LOAD_X, dataSend), d);
                } else {
                    mpi.send(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.LOAD_W, dataSend), d);
                }
            }
        } else { // parallel
            ParallelLoading[] pl = new ParallelLoading[domain.nDoms];
            for (int d = 0; d < domain.nDoms; d++) {
                pl[d] = new ParallelLoading(d, index);
                pl[d].start();
            }
            try {
                for (int d = 0; d < domain.nDoms; d++) {
                    pl[d].join();
                }
            } catch (java.lang.InterruptedException e) {
                System.out.println(e);
            }
        }

        mpi.waitForAll(Tag.GMRES2MASTER);
    }

    double norm(char s) throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.NORM, s));
        double sum = mpi.receiveAllDoubleSum(Tag.GMRES2MASTER);
        return Math.sqrt(sum);
    }

    double[] timesVectorByScalar(double[] a, double b) {
        double[] c = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            c[i] = a[i] * b;
        }
        return c;
    }

    double[] rotmat(double a, double b) {
        // Compute the Givens rotation matrix parameters for a and b.
        double c, s, temp;
        if (b == 0) {
            c = 1;
            s = 0;
        } else if (Math.abs(b) > Math.abs(a)) {
            temp = a / b;
            s = 1 / Math.sqrt(1 + temp * temp);
            c = temp * s;
        } else {
            temp = b / a;
            c = 1 / Math.sqrt(1 + temp * temp);
            s = temp * c;
        }
        return new double[]{c, s};
    }

    // Gaussian elimination with partial pivoting
    double[] lsolve(double[][] A, double[] b, int N) {
        // int N  = b.length;
        for (int p = 0; p < N; p++) {
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        // back substitution
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

    class ParallelSaving extends Thread {

        int d;

        ParallelSaving(int d) {
            this.d = d;
        }

        public void run() {
            try {
                int[] mapL2G = domain.getSubdomain(d).mapL2G;
                LiteElement[] dataRcv = (LiteElement[]) mpi.receive(d, Tag.GMRES2MASTER).getData();
                for (LiteElement dataRcv1 : dataRcv) {
                    //liteElems[mapL2G[dataRcv1.index]] = new LiteElement(dataRcv1.index, dataRcv1.y);
                    liteElems[mapL2G[dataRcv1.index]] = dataRcv1;
                }
            } catch (MPIException ex) {
                Logger.getLogger(ParallelGmresMaster.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    class ParallelLoading extends Thread {

        int d;
        int index;

        ParallelLoading(int d, int index) {
            this.d = d;
            this.index = index;
        }

        public void run() {
            try {
                Domain.Subdomain subdom = domain.getSubdomain(d);
                int[] mapG2L = subdom.mapG2L;
                int[] load = subdom.load1;
                LiteElement[] dataSend = new LiteElement[load.length];
                for (int i = 0; i < load.length; i++) {
                    //dataSend[i] = new LiteElement(mapG2L[load[i]], liteElems[load[i]].y);
                    dataSend[i] = liteElems[load[i]];
                    dataSend[i].index = mapG2L[load[i]];
                }
                if (index == -1) {
                    mpi.send(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.LOAD_X, dataSend), d);
                } else {
                    mpi.send(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.LOAD_W, dataSend), d);
                }
            } catch (MPIException ex) {
                Logger.getLogger(ParallelGmresMaster.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void ping() throws MPIException {
        Watch w = new Watch();
        w.start();
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, ParallelTags.PING));
        mpi.waitForAll(Tag.GMRES2MASTER);
        w.stop();
    }
}
