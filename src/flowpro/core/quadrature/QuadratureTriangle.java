package flowpro.core.quadrature;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * Quadrature1DLine for integrating functions over triangles or rectangles.
 *
 * @author ales
 */
public class QuadratureTriangle extends Quadrature {

    public QuadratureTriangle(String fileName) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            dimension = 2;
            nPoints = Integer.parseInt(reader.readLine());
            weights = new double[nPoints];

            double[][] barCoords = new double[nPoints][3];
            for (int i = 0; i < nPoints; ++i) {
                String line = reader.readLine();
                StringTokenizer tokenizer = new StringTokenizer(line);
                for (int j = 0; j < 3; ++j) {
                    barCoords[i][j] = Double.parseDouble(tokenizer.nextToken());
                }
            }
            // Cartesian coordinates of the reference triangle
            double[] x = {0, 1, 0};
            double[] y = {0, 0, 1};
            coords = setCoords(barCoords, x, y);

            for (int i = 0; i < nPoints; ++i) {
                weights[i] = Double.parseDouble(reader.readLine());
            }
        } catch (NullPointerException | NumberFormatException |
                ArrayIndexOutOfBoundsException ex) {
            throw new IOException("reading quadrature points and weights failed: file "
                    + fileName + " has a wrong format");
        } catch (IOException ex) {
            throw new IOException("reading quadrature points and weights failed: "
                    + ex.getMessage());
        }
    }
    
    public double[][] setCoords(double[][] barCoords, double[] x, double[] y) {
        int nPoints = barCoords.length;
        double[][] coords = new double[nPoints][2];

        for (int i = 0; i < nPoints; ++i) {
            for (int j = 0; j < 3; ++j) {
                coords[i][0] += x[j] * barCoords[i][j];
                coords[i][1] += y[j] * barCoords[i][j];
            }
        }

        return coords;
    }

    public double[][] getCoords() { 
        return coords;
    }
}
