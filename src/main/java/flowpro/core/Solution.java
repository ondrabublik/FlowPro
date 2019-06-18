package flowpro.core;

import java.io.*;
import flowpro.core.parallel.Domain;

/**
 *
 * @author ales
 */
public class Solution implements Serializable {

    public final int nElems;
//    public final int nPoints;

    public final double[][] W;
    public final double[][] avgW;
    public double[][] vertices = null;
    public double[][] meshVelocity = null;

    public Solution(double[][] W, double[][] avgW) {
        nElems = W.length;

        this.W = W;
        this.avgW = avgW;
    }

    public Solution(double[][] W, double[][] avgW, double[][] PXY) {
        nElems = W.length;

        this.W = W;
        this.avgW = avgW;
        this.vertices = PXY;
    }

    public Solution(double[][] W, double[][] avgW, double[][] PXY, double[][] U) {
        nElems = W.length;

        this.W = W;
        this.avgW = avgW;
        this.vertices = PXY;
        this.meshVelocity = U;
    }

    public Solution(Solution[] sols, Domain domain) {
        nElems = domain.nElems;
        int nPoints = domain.nPoints;
        W = new double[nElems][];
        avgW = new double[nElems][];
        int dim = sols[0].vertices[0].length;
        vertices = new double[nPoints][dim];
        if (sols[0].meshVelocity != null) {
            meshVelocity = new double[nPoints][dim];
        }

        for (int d = 0; d < sols.length; ++d) {
            Domain.Subdomain subdom = domain.getSubdomain(d);
            int[][] TP = subdom.TP;
            int[] mapL2G = subdom.mapL2G;
            for (int i : subdom.interior) {
                int pos = mapL2G[i];
                W[pos] = sols[d].W[i];
                avgW[pos] = sols[d].avgW[i];
                for (int j = 0; j < TP[i].length; j++) {
                    for (int k = 0; k < dim; k++) {
                        vertices[TP[i][j]][k] = sols[d].vertices[TP[i][j]][k];
                    }
                }
                if (sols[d].meshVelocity != null) {
                    for (int j = 0; j < TP[i].length; j++) {
                        for (int k = 0; k < dim; k++) {
                            meshVelocity[TP[i][j]][k] = sols[d].meshVelocity[TP[i][j]][k];
                        }
                    }
                }
            }
        }
    }

    public void saveSolution(String fileName) {
        try {
            FileWriter fw = new FileWriter(fileName);
            try (BufferedWriter out = new BufferedWriter(fw)) {
                appendDoubleMatrix(out, avgW, "avgW");
//                appendDoubleMatrix(out, detailW, "detailW");
                appendDoubleMatrix(out, W, "W");
                appendDoubleMatrix(out, vertices, "PXY");
            }
        } catch (IOException e) {
            System.out.println(e);
        }
    }

    private void appendDoubleMatrix(BufferedWriter out, double[][] M, String str) throws IOException {
        //if(M != null){
        //   String radka = str + " " + M.length;

        String radka = str + " " + M.length;
        out.write(radka);
        out.newLine();
        for (double[] M1 : M) {
            radka = M1.length + " ";
            for (int j = 0; j < M1.length; j++) {
                radka = radka + " " + String.valueOf(M1[j]);
            }
            out.write(radka);
            out.newLine();
            for (int i = 0; i < M.length; i++) {
                radka = M[i].length + " ";
                for (int j = 0; j < M[i].length; j++) {
                    radka = radka + " " + String.valueOf(M[i][j]);
                }
                out.write(radka);
                out.newLine();
            }
        }
    }
}
