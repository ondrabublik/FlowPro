package flowpro.core.quadrature;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * QuadratureLine for integrating scalar single valued functions.
 *
 * @author ales
 */
public class QuadratureLine extends Quadrature {
    
    public QuadratureLine(String fileName) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            dimension = 1;
            nPoints = Integer.parseInt(reader.readLine());
            coords = new double[nPoints][dimension];
            weights = new double[nPoints];

            for (int i = 0; i < nPoints; ++i) {
                coords[i][0] = Double.parseDouble(reader.readLine());
            }

            for (int i = 0; i < nPoints; ++i) {
                weights[i] = Double.parseDouble(reader.readLine());
            }
        } catch (NullPointerException | NumberFormatException ex) {
            throw new IOException("reading quadrature points and weights failed: file "
                    + fileName + " has a wrong format");
        } catch (IOException ex) {
            throw new IOException("reading quadrature points and weights failed: " + ex.getMessage());
        }
    }
    
    public double[][] getCoords(){
        return coords;
    }
}
