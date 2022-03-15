package flowpro.core.quadrature;

import java.io.IOException;
import java.io.Serializable;

/**
 *
 * @author obublik
 */
public class QuadratureCentral implements Serializable {

    public Quadrature[] qPoint;
    public Quadrature[] qLine;
    public Quadrature[] qTriangle;
    public Quadrature[] qSquare;
    public Quadrature[] qTetra;
    public Quadrature[] qHexa;
    public Quadrature[] qPrism;
    public Quadrature[] q25DSquare;

    public QuadratureCentral() throws IOException {
        int orderMax = 10;

        // load Gaussian quadrature
        qPoint = new Quadrature[orderMax + 1];
        qLine = new Quadrature[orderMax + 1];
        qTriangle = new Quadrature[orderMax + 1];
        qSquare = new Quadrature[orderMax + 1];
        qTetra = new Quadrature[orderMax + 1];
        qHexa = new Quadrature[orderMax + 1];
        qPrism = new Quadrature[orderMax + 1];
        q25DSquare = new Quadrature[orderMax + 1];

        for (int order = 1; order <= orderMax; order++) {
            qPoint[order] = new QuadraturePoint();
            try {
                String quadratureFileName = "gaussRules/gauss/" + order + ".txt";
                qLine[order] = new QuadratureLine(quadratureFileName);
            } catch (Exception e) {
            }
            try {
                int ord = 2 * (order - 1);
                if(ord == 0){
                    ord = 1;
                }
                String quadratureFileName = "gaussRules/gaussTri/" + ord + ".txt";
                qTriangle[order] = new QuadratureTriangle(quadratureFileName);
            } catch (Exception e) {
            }
            try {
                String quadratureFileName = "gaussRules/gauss/" + order + ".txt";
                qSquare[order] = new QuadratureSquare(quadratureFileName);
            } catch (Exception e) {
            }
            try {
                int ord = 2 * (order - 1);
                if(ord == 0){
                    ord = 1;
                }
                String quadratureFileName = "gaussRules/gaussTetra/" + ord + ".txt";
                qTetra[order] = new QuadratureTetra(quadratureFileName);
            } catch (Exception e) {
            }
            try {
                String quadratureFileName = "gaussRules/gauss/" + order + ".txt";
                qHexa[order] = new QuadratureHexa(quadratureFileName);
            } catch (Exception e) {
            }
            try {
                qPrism[order] = new QuadraturePrism(qTriangle[order], qLine[order]);
            } catch (Exception e) {
            }
            try {
                q25DSquare[order] = new QuadratureSquare25D(order);
            } catch (Exception e) {
                System.out.println(e);
            }
        }
    }
}
