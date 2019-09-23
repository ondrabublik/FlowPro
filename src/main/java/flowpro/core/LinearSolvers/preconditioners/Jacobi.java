/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.core.LinearSolvers.SparseMatrix;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import flowpro.core.element.Implicit;

/**
 *
 * @author obublik
 */
class Jacobi extends Preconditioner {

    SparseMatrix A;
    int n, nThreads;

    Jacobi(Parameters par) {
        nThreads = par.nThreads;
    }

    @Override
    public void setMatrix(SparseMatrix A) {
        this.A = A;
        n = A.getDofs();
    }

    @Override
    public void factor() {

    }

    @Override
    public void apply(double[] x, double[] b) {

        Element[] elems = A.getElems();

        // vlastni vypocet, parallelni beh
        Pthread[] parallel = new Pthread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new Pthread(v, nThreads, x, b, elems);
            parallel[v].start();
        }
        try {
            for (int v = 0; v < nThreads; v++) {
                parallel[v].join();
            }
        } catch (java.lang.InterruptedException e) {
            System.out.println(e);
        }
    }

    class Pthread extends Thread {

        int nStart, nThreads, par, nt;
        double[] x, b;
        Element[] elems;

        Pthread(int nStart, int nThreads, double[] x, double[] b, Element[] elems) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.x = x;
            this.b = b;
            this.elems = elems;
            nt = elems.length;
        }

        @Override
        public void run() {
            for (int i = nStart; i < nt; i = i + nThreads) {
                Element elem = elems[i];
                if (elem.insideComputeDomain) {
                    double[][] Adiag = ((Implicit)elem.ti).ADiag;
                    int[] glob = elem.gi_U;
                    for (int j = 0; j < glob.length; j++) {
                        x[glob[j]] = 1 / Adiag[j][j] * b[glob[j]];
                    }
                }
            }
        }
    }
}
