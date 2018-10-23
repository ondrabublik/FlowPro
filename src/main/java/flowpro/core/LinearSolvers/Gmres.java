/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package flowpro.core.LinearSolvers;

import flowpro.core.LinearSolvers.preconditioners.Preconditioner;

/**
 *
 * @author obublik
 */
public class Gmres extends LinearSolver{

    int n, m, iterationLimit, nThreads;
    double tol;
    double[][] V, H;
    double[] cs, sn, e1, w, r, aux;

    Gmres(SparseMatrix A, Preconditioner M, int m, int iterationLimit, double tol, int nThreads) {
        this.A = A;
        A.buildCRSformat();
        this.M = M;
        this.n = A.getDofs();
        this.m = m;
        this.iterationLimit = iterationLimit;
        this.tol = tol;
        this.nThreads = nThreads;

        // initialize workspace
        V = new double[m + 1][n];
        H = new double[m + 1][m];
        cs = new double[m];
        sn = new double[m];
        e1 = new double[m + 1];
        w = new double[n];
        r = new double[n];
        aux = new double[n];
    }

    public boolean solve(double[] x, double[] b) {
        double norm_r, temp;
        double[] s, y;

        double bnrm2 = 0;
        for (int i = 0; i < n; i++) {
            bnrm2 += b[i] * b[i];
        }
        bnrm2 = Math.sqrt(bnrm2);

        if (bnrm2 == 0) {
            bnrm2 = 1;
        }

        A.SubstrMult(aux, b, x, nThreads);
        M.apply(r, aux);
        double error = norm(r) / bnrm2;

        if (error < tol) {
            return true;
        }
        e1[0] = 1;

        for (int iter = 0; iter < iterationLimit; iter++) {          // begin iteration
            norm_r = norm(r);
            for (int j = 0; j < n; j++) {
                V[0][j] = r[j] / norm_r;
            }
            s = vectorScalarProduct(e1, norm_r);
            for (int i = 0; i < m; i++) {                        // construct orthonormal
                A.Mult(aux, V[i], nThreads);     // basis using Gram-Schmidt
                M.apply(w, aux);
                for (int k = 0; k <= i; k++) {
                    H[k][i] = scalarProduct(w, V[k]);
                    for (int j = 0; j < n; j++) {
                        w[j] = w[j] - V[k][j] * H[k][i];
                    }
                }
                H[i + 1][i] = norm(w);
                for (int j = 0; j < n; j++) {
                    V[i + 1][j] = w[j] / H[i + 1][i];
                }
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
                error = Math.abs(s[i + 1]) / bnrm2;

                if (error <= tol) {                        // update approximation
                    y = lsolve(H, s, i + 1);                 // and exit
                    updateSolution(y, x, i, nThreads);
                    break;
                }
            }

            if (error <= tol) {
                break;
            }

            y = lsolve(H, s, m);
            updateSolution(y, x, m - 1, nThreads);
            A.SubstrMult(aux, b, x, nThreads);                     // compute residual
            M.apply(r, aux);
            s[m] = norm(r);
            error = s[m] / bnrm2;                     // check convergence
            if (error <= tol) {
                break;
            }
        }
        return error <= tol;
    }

    void updateSolution(double[] y, double[] x, int i, int nThreads) {
        // vlastni vypocet, parallelni beh
        UpdateSolutionThread[] parallel = new UpdateSolutionThread[nThreads];
        for (int v = 0; v < nThreads; v++) {
            parallel[v] = new UpdateSolutionThread(v, nThreads, y, x, i);
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

    class UpdateSolutionThread extends Thread {

        int nStart, nThreads, i, n;
        double[] y, x;

        UpdateSolutionThread(int nStart, int nThreads, double[] y, double[] x, int i) {
            this.nStart = nStart;
            this.nThreads = nThreads;
            this.y = y;
            this.x = x;
            this.i = i;
            n = x.length;
        }

        @Override
        public void run() {
            // paralelni sestavovani a plneni matic
            for (int j = nStart; j < n; j = j + nThreads) {
                for (int k = 0; k <= i; k++) {
                    x[j] = x[j] + V[k][j] * y[k];
                }
            }
        }
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

    double[] vectorScalarProduct(double[] a, double b) {
        double[] c = new double[n];
        for (int i = 0; i < a.length; i++) {
            c[i] = a[i] * b;
        }
        return c;
    }

    double norm(double[] a) {
        double n = 0;
        for (int i = 0; i < a.length; i++) {
            n = n + a[i] * a[i];
        }
        n = Math.sqrt(n);

        return n;
    }

    double scalarProduct(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            s = s + a[i] * b[i];
        }
        return s;
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
}
