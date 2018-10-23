//package flowpro.core.LinearSolvers;
//
//import flowpro.core.Mesh.Element;
//
///**
// *
// * @author obublik
// */
//public class BiCgStab extends LinearSolver2 {
//
//    int n;
//    int iterationLimit;
//    int nThreads;
//    double tol;
//
//    double[] r, s, r0, p, Ap, As;
//    double alfa, omega, beta, rr0, rr0n;
//
//    public BiCgStab(Element[] elems, int n, int iterationLimit, double tol, int nThreads) {
//        this.elems = elems;
//        this.n = n;
//        this.iterationLimit = iterationLimit;
//        this.tol = tol;
//        this.nThreads = nThreads;
//
//        r = new double[n];
//        s = new double[n];
//        Ap = new double[n];
//        As = new double[n];
//    }
//
//    @Override
//    public boolean solve(double[] x) {
//
//        ComputeResiduum(x, r, 1, nThreads);
//        r0 = copy(r);
//        p = copy(r);
//        rr0 = scalarProduct(r, r0);
//
//        for (int iter = 0; iter < iterationLimit; iter++) {
//            ComputeResiduum(p, Ap, 0, nThreads);
//            alfa = rr0 / scalarProduct(Ap, r0);
//            for (int i = 0; i < n; i++) {
//                s[i] = r[i] - alfa * Ap[i];
//            }
//            ComputeResiduum(s, As, 0, nThreads);
//            omega = scalarProduct(As, s) / scalarProduct(As, As);
//            for (int i = 0; i < n; i++) {
//                x[i] = x[i] + alfa * p[i] + omega * s[i];
//                r[i] = s[i] - omega * As[i];
//            }
//            rr0n = scalarProduct(r, r0);
//            beta = rr0n / rr0 * alfa / omega;
//            rr0 = rr0n;
//            for (int i = 0; i < n; i++) {
//                p[i] = r[i] + beta*(p[i] - omega*Ap[i]);
//            }
//            
//            double error = 0;
//            for (int i = 0; i < n; i++) {
//                error += r[i]*r[i];
//            }
//            
//            if(Math.sqrt(error) < tol){
//                return true;
//            }
//        }
//
//        return false;
//    }
//    
//    void ComputeResiduum(double[] x, double[] r, int par, int nThreads) {
//        // vlastni vypocet, parallelni beh
//        BiCgStabThread[] parallel = new BiCgStabThread[nThreads];
//        for (int v = 0; v < nThreads; v++) {
//            parallel[v] = new BiCgStabThread(v, nThreads, x, r, par);
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
//    class BiCgStabThread extends Thread {
//
//        int nStart, nThreads, par, nt;
//        double[] x, r;
//        double residuumJacobi;
//
//        BiCgStabThread(int nStart, int nThreads, double[] x, double[] r, int par) {
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
