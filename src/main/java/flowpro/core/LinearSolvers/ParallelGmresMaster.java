/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.Parameters;
import flowpro.core.parallel.Domain;
import flowpro.core.parallel.LiteElement;
import flowpro.core.parallel.Tag;
import litempi.*;

/**
 *
 * @author obublik
 */
public class ParallelGmresMaster {

    int iterationLimit, nThreads;
    double tol;
    MPIMaster mpi;
    Domain domain;
    LiteElement[] liteElems;

    public ParallelGmresMaster(Parameters par, MPIMaster mpi, Domain domain) {
        this.mpi = mpi;
        this.domain = domain;
        iterationLimit = 300;
        tol = par.iterativeSolverTol;
        nThreads = par.nThreads;
        liteElems = new LiteElement[domain.nElems];
    }

    public boolean solve() throws MPIException {
        boolean converged = false;
        for (int i = 0; i < iterationLimit; i++) {
            multiplyAndUpdate();
            dataExchange();
            double error = norm();
            if (error < tol) {
                converged = true;
                break;
            }
        }
        return converged;
    }

    void multiplyAndUpdate() throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, 0));
        mpi.waitForAll(Tag.GMRES2MASTER);
    }

    void dataExchange() throws MPIException {
        // downloading data from the slave nodes into the central structure
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE,3));
        for (int d = 0; d < domain.nDoms; ++d) {
            int[] mapL2G = domain.getSubdomain(d).mapL2G;
            LiteElement[] dataRcv = (LiteElement[]) mpi.receive(d, Tag.GMRES2MASTER).getData();
            for (LiteElement dataRcv1 : dataRcv) {
                liteElems[mapL2G[dataRcv1.index]] = new LiteElement(dataRcv1.index, dataRcv1.y);
            }
        }

        // uploading data from the central structure onto the slave nodes
        for (int d = 0; d < domain.nDoms; ++d) {
            Domain.Subdomain subdom = domain.getSubdomain(d);
            int[] mapG2L = subdom.mapG2L;
            int[] load = subdom.load1;
            LiteElement[] dataSend = new LiteElement[load.length];
            for (int i = 0; i < load.length; i++) {
                dataSend[i] = new LiteElement(mapG2L[load[i]], liteElems[load[i]].y);
            }
            mpi.send(new MPIMessage(Tag.GMRES2SLAVE, 4, dataSend), d);
        }
        mpi.waitForAll(Tag.GMRES2MASTER);
    }

    double norm() throws MPIException {
        mpi.sendAll(new MPIMessage(Tag.GMRES2SLAVE, 2));
        return mpi.receiveAllDoubleSum(Tag.GMRES2MASTER);
    }
}
