//package flowpro.core.LinearSolvers;
//
//import flowpro.core.Mesh.*;
//import no.uib.cipr.matrix.DenseVector;
//import no.uib.cipr.matrix.Vector;
//import no.uib.cipr.matrix.sparse.BiCGstab;
//import no.uib.cipr.matrix.sparse.CompRowMatrix;
//import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
//import no.uib.cipr.matrix.sparse.GMRES;
//import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
//
///**
// *
// * @author obublik
// */
//public class MTJsolver extends LinearSolver {
//
//    int n;
//    int iterationLimit;
//    int nThreads;
//    double tol;
//    //GMRES solver;
//    BiCGstab solver;
//
//    public MTJsolver(Element[] elems, int n, int iterationLimit, double tol, int nThreads) {
//        this.elems = elems;
//        this.n = n;
//        this.iterationLimit = iterationLimit;
//        this.tol = tol;
//        this.nThreads = nThreads;
//        System.out.println("Hello to MTJ solver!");
//        //AMG prec = new AMG();
//        //ILU prec = new ILU(A);
//        //prec.setMatrix(A);
//        //solver.setPreconditioner(prec);   
//    }
//
//    @Override
//    public boolean solve(double[] x) {
//        try {
//            FlexCompRowMatrix A = new FlexCompRowMatrix(n, n);
//            Vector b = new DenseVector(n);
//            for (int el = 0; el < elems.length; el++) {
//                for (int k = 0; k < elems[el].nFaces; k++) {
//                    if (elems[el].TT[k] > -1) {
//                        for (int i = 0; i < elems[el].gi_U.length; i++) {
//                            for (int j = 0; j < elems[elems[el].TT[k]].gi_U.length; j++) {
//                                A.add(elems[el].gi_U[i],elems[elems[el].TT[k]].gi_U[j],  elems[el].ANeighs[k].MR[j][i]);
//                            }
//                        }
//                    }
//                }
//                for (int i = 0; i < elems[el].gi_U.length; i++) {
//                    b.add(elems[el].gi_U[i], elems[el].RHS_loc[i]);
//                    for (int j = 0; j < elems[el].gi_U.length; j++) {
//                        A.add(elems[el].gi_U[i], elems[el].gi_U[j], elems[el].ADiag[i][j]);
//                    }
//                }
//            }
//            CompRowMatrix C = new CompRowMatrix(A);
//            Vector y = new DenseVector(n);
//            //solver = new GMRES(b,50);
//            solver = new BiCGstab(b);
////            ILU prec = new ILU(C);
////            prec.setMatrix(C);
////            solver.setPreconditioner(prec); 
//            solver.solve(C, b, y);
//            for (int i = 0; i < x.length; i++) {
//                x[i] = y.get(i);
//            }
//
//            return true;
//        } catch (IterativeSolverNotConvergedException ex) {
//            return false;
//        }
//    }
//}
