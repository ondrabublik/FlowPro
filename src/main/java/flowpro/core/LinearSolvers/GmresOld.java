//package flowpro.core.LinearSolvers;
//
//import flowpro.core.Mesh.Element;
//
///**
// *
// * @author ales
// */
//public class Gmres extends LinearSolver {
//
//    int n;
//    int m;
//    int iterationLimit;
//    int nThreads;
//    double tol;
//    double[][] V, H;
//    double[] cs, sn, e1, w, r;
//
//    public Gmres(Element[] elems, int n, int m, int iterationLimit, double tol, int nThreads) {
//        this.elems = elems;
//        this.n = n;
//        this.m = m;
//        this.iterationLimit = iterationLimit;
//        this.tol = tol;
//        this.nThreads = nThreads;
//
//        // initialize workspace
//        V = new double[m + 1][n];
//        H = new double[m + 1][m];
//        cs = new double[m];
//        sn = new double[m];
//        e1 = new double[m + 1];
//        w = new double[n];
//        r = new double[n];
//    }
//
//    @Override
//    public boolean solve(double[] x) {
//        double norm_r, temp;
//        double[] s, y;
//
//        double bnrm2 = 0;
//        for (Element elem : elems) {
//            bnrm2 = bnrm2 + elem.sqr();
//        }
//        bnrm2 = Math.sqrt(bnrm2);
//
//        if (bnrm2 == 0) {
//            bnrm2 = 1;
//        }
//
//        ComputeResiduum(x, r, 1, nThreads);
//        double error = norm(r) / bnrm2;
//        if (error < tol) {
//            return true;
//        }
//        e1[0] = 1;
//
//        for (int iter = 0; iter < iterationLimit; iter++) {          // begin iteration
//            //ComputeResiduum(x, r, 1, nThreads);
//            norm_r = norm(r);
//            for (int j = 0; j < n; j++) {
//                V[0][j] = r[j] / norm_r;
//            }
//            s = vectorScalarProduct(e1, norm_r);
//            for (int i = 0; i < m; i++) {                        // construct orthonormal
//                ComputeResiduum(V[i], w, 0, nThreads);     // basis using Gram-Schmidt
//                for (int k = 0; k <= i; k++) {
//                    H[k][i] = scalarProduct(w, V[k]);
//                    for (int j = 0; j < n; j++) {
//                        w[j] = w[j] - V[k][j] * H[k][i];
//                    }
//                }
//                H[i + 1][i] = norm(w);
//                for (int j = 0; j < n; j++) {
//                    V[i + 1][j] = w[j] / H[i + 1][i];
//                }
//                for (int k = 0; k <= i - 1; k++) {                              // apply Givens rotation
//                    temp = cs[k] * H[k][i] + sn[k] * H[k + 1][i];
//                    H[k + 1][i] = -sn[k] * H[k][i] + cs[k] * H[k + 1][i];
//                    H[k][i] = temp;
//                }
//                double[] rot = rotmat(H[i][i], H[i + 1][i]); // form i-th rotation matrix
//                cs[i] = rot[0];
//                sn[i] = rot[1];
//                temp = cs[i] * s[i];                            // approximate residual norm
//                s[i + 1] = -sn[i] * s[i];
//                s[i] = temp;
//                H[i][i] = cs[i] * H[i][i] + sn[i] * H[i + 1][i];
//                H[i + 1][i] = 0;
//                error = Math.abs(s[i + 1]) / bnrm2;
//                if (error <= tol) {                        // update approximation
//                    y = lsolve(H, s, i + 1);                 // and exit
//                    updateSolution(V, y, x, i, nThreads);
//                    break;
//                }
//            }
//            if (error <= tol) {
//                break;
//            }
//
//            y = lsolve(H, s, m);
//            updateSolution(V, y, x, m - 1, nThreads);
//            ComputeResiduum(x, r, 1, nThreads);                      // compute residual
//            s[m] = norm(r);
//            error = s[m] / bnrm2;                     // check convergence
//            if (error <= tol) {
//                break;
//            }
//        }
//        return error <= tol;
//    }
//
//    // presunout do tridy MAT !!!!!!!
//    void updateSolution(double[][] V, double[] y, double[] x, int i, int nThreads) {
//        // vlastni vypocet, parallelni beh
//        UpdateSolutionThread[] parallel = new UpdateSolutionThread[nThreads];
//        for (int v = 0; v < nThreads; v++) {
//            parallel[v] = new UpdateSolutionThread(v, nThreads, V, y, x, i);
//            parallel[v].start();
//        }
//        try {
//            for (int v = 0; v < nThreads; v++) {
//                parallel[v].join();
//            }
//        } catch (java.lang.InterruptedException e) {
//            System.out.println(e);
//        }
//    }
//
//    double[] rotmat(double a, double b) {
//        // Compute the Givens rotation matrix parameters for a and b.
//        double c, s, temp;
//        if (b == 0) {
//            c = 1;
//            s = 0;
//        } else if (Math.abs(b) > Math.abs(a)) {
//            temp = a / b;
//            s = 1 / Math.sqrt(1 + temp * temp);
//            c = temp * s;
//        } else {
//            temp = b / a;
//            c = 1 / Math.sqrt(1 + temp * temp);
//            s = temp * c;
//        }
//        return new double[]{c, s};
//    }
//
//    double[] vectorScalarProduct(double[] a, double b) {
//        double[] c = new double[n];
//        for (int i = 0; i < a.length; i++) {
//            c[i] = a[i] * b;
//        }
//        return c;
//    }
//
//    double norm(double[] a) {
//        double n = 0;
//        for (int i = 0; i < a.length; i++) {
//            n = n + a[i] * a[i];
//        }
//        n = Math.sqrt(n);
//
//        return n;
//    }
//
//    // Gaussian elimination with partial pivoting
//    double[] lsolve(double[][] A, double[] b, int N) {
//        // int N  = b.length;
//        for (int p = 0; p < N; p++) {
//            // find pivot row and swap
//            int max = p;
//            for (int i = p + 1; i < N; i++) {
//                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
//                    max = i;
//                }
//            }
//            double[] temp = A[p];
//            A[p] = A[max];
//            A[max] = temp;
//            double t = b[p];
//            b[p] = b[max];
//            b[max] = t;
//
//            // pivot within A and b
//            for (int i = p + 1; i < N; i++) {
//                double alpha = A[i][p] / A[p][p];
//                b[i] -= alpha * b[p];
//                for (int j = p; j < N; j++) {
//                    A[i][j] -= alpha * A[p][j];
//                }
//            }
//        }
//
//        // back substitution
//        double[] x = new double[N];
//        for (int i = N - 1; i >= 0; i--) {
//            double sum = 0.0;
//            for (int j = i + 1; j < N; j++) {
//                sum += A[i][j] * x[j];
//            }
//            x[i] = (b[i] - sum) / A[i][i];
//        }
//        return x;
//    }
//    
//    void ComputeResiduum(double[] x, double[] r, int par, int nThreads) {
//        // vlastni vypocet, parallelni beh
//        GmresThread[] parallel = new GmresThread[nThreads];
//        for (int v = 0; v < nThreads; v++) {
//            parallel[v] = new GmresThread(v, nThreads, x, r, par);
//            parallel[v].start();
//        }
//        try {
//            for (int v = 0; v < nThreads; v++) {
//                parallel[v].join();
//            }
//        } catch (java.lang.InterruptedException e) {
//            System.out.println(e);
//        }
//    }
//
//    class GmresThread extends Thread {
//
//        int nStart, nThreads, par, nt;
//        double[] x, r;
//        double residuumJacobi;
//
//        GmresThread(int nStart, int nThreads, double[] x, double[] r, int par) {
//            this.nStart = nStart;
//            this.nThreads = nThreads;
//            this.x = x;
//            this.r = r;
//            this.par = par;
//            nt = elems.length;
//        }
//
//        @Override
//        public void run() {
//            for (int i = nStart; i < nt; i = i + nThreads) {
//                if (elems[i].insideComputeDomain) {
//                    elems[i].residuumGmres(x, r, par);
//                }
//            }
//        }
//    }
//}
//
//class UpdateSolutionThread extends Thread {
//
//    int nStart, nThreads, i, n;
//    double[][] V;
//    double[] y, x;
//
//    UpdateSolutionThread(int nStart, int nThreads, double[][] V, double[] y, double[] x, int i) {
//        this.nStart = nStart;
//        this.nThreads = nThreads;
//        this.V = V;
//        this.y = y;
//        this.x = x;
//        this.i = i;
//        n = x.length;
//    }
//
//    @Override
//    public void run() {
//        // paralelni sestavovani a plneni matic
//        for (int j = nStart; j < n; j = j + nThreads) {
//            for (int k = 0; k <= i; k++) {
//                x[j] = x[j] + V[k][j] * y[k];
//            }
//        }
//    }
//}
