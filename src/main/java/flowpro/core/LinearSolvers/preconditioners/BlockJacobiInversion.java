/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers.preconditioners;

import flowpro.api.Mat;
import flowpro.core.LinearSolvers.SparseMatrix;
import flowpro.core.Mesh.Element;
import flowpro.core.Parameters;

/**
 *
 * @author obublik
 */
class BlockJacobiInversion extends Preconditioner {

    SparseMatrix A;
    int n, nThreads;
    double[][][] diagonalInverse;
    Element[] elems;

    BlockJacobiInversion(Parameters par) {
        nThreads = par.nThreads;
    }

    @Override
    public void setMatrix(SparseMatrix A) {
        this.A = A;
        n = A.getDofs();
        diagonalInverse = new double[n][][];
        elems = A.getElems();
    }

    @Override
    public void factor() {
        // vlastni vypocet, parallelni beh
        FactorThread[] parallel = new FactorThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new FactorThread(v, nThreads);
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

    class FactorThread extends Thread {

        int nStart, nThreads, par;
        double[] x, b;

        FactorThread(int nStart, int nThreads) {
            this.nStart = nStart;
            this.nThreads = nThreads;
        }

        @Override
        public void run() {
            for (int i = nStart; i < elems.length; i = i + nThreads) {
                if (elems[i].insideComputeDomain) {
                    diagonalInverse[i] = Mat.invert(elems[i].ADiag);
                }
            }
        }
    }

    @Override
    public void apply(double[] x, double[] b) {

        // vlastni vypocet, parallelni beh
        ApplyThread[] parallel = new ApplyThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new ApplyThread(v, nThreads, x, b);
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

    class ApplyThread extends Thread {

        int nStart, nThreads, par;
        double[] x, b;

        ApplyThread(int nStart, int nThreads, double[] x, double[] b) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.x = x;
            this.b = b;
        }

        @Override
        public void run() {
            for (int i = nStart; i < elems.length; i = i + nThreads) {
                Element elem = elems[i];
                if (elem.insideComputeDomain) {
                    int[] glob = elem.gi_U;
                    double[][] M = diagonalInverse[i];
                    for (int j = 0; j < glob.length; j++) {
                        x[glob[j]] = 0;
                        for (int k = 0; k < glob.length; k++) {
                            x[glob[j]] += M[k][j] * b[glob[k]];
                        }
                    }
                }
            }
        }
    }
}
