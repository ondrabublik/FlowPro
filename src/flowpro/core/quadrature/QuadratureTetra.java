package flowpro.core.quadrature;

import flowpro.api.Mat;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * Quadrature1DLine for integrating functions over triangles or rectangles.
 *
 * @author ales
 */
public class QuadratureTetra extends Quadrature {

    public QuadratureTetra(String fileName) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            dimension = 3;
            nPoints = Integer.parseInt(reader.readLine());
            weights = new double[nPoints];

            coords = new double[nPoints][dimension];
            for (int i = 0; i < nPoints; ++i) {
                String line = reader.readLine();
                StringTokenizer tokenizer = new StringTokenizer(line);
                for (int j = 0; j < 3; ++j) {
                    coords[i][j] = Double.parseDouble(tokenizer.nextToken());
                }
            }
            
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

    public double[][] getCoords() { 
        return coords;
    }
}