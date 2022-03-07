
package flowpro.core.basis;

import java.io.IOException;

/**
 * @author obublik
**/

public class basis3D_25D extends Basis {
 
    private int order;

    public basis3D_25D(int order) throws IOException {
        this.order = order;
        int[] nBas = new int[]{1, 3, 9};
        nBasis = nBas[order - 1];
        
        
    }

    @Override
    public void calculateCoefficients() {
        basisType = "taylor";
    }

    // tato funkce vraci hodnotu m-te bazove funkce v bode Xi
    @Override
    public double basisFun(int m, double[] Xi) {
        double val = 0;
        switch(m) {
            case 0:
                val = 1;
                break;
            case 1:
                val = Xi[0];
                break;
            case 2:
                val = Xi[1];
                break;
            case 3:
                val = Math.sin(2*Math.PI*Xi[2]);
                break;
            case 4:
                val = Xi[0]*Math.sin(2*Math.PI*Xi[2]);;
                break;
            case 5:
                val = Xi[1]*Math.sin(2*Math.PI*Xi[2]);;
                break;
            case 6:
                val = Math.cos(2*Math.PI*Xi[2]);
                break;
            case 7:
                val = Xi[0]*Math.cos(2*Math.PI*Xi[2]);;
                break;
            case 8:
                val = Xi[1]*Math.cos(2*Math.PI*Xi[2]);;
                break;
        }
        return val;
    }

    // tato funkce vraci hodnotu x-ove derivace m-te bazove funkce v bode x, y
    @Override
    public double derBasis(int m, double[] Xi, int dim) {
        double val = 0;
        double c = 2*Math.PI;
        switch (dim) {
            case 0:
                switch(m) {
                    case 0:
                        val = 0;
                        break;
                    case 1:
                        val = 1;
                        break;
                    case 2:
                        val = 0;
                        break;
                    case 3:
                        val = 0;
                        break;
                    case 4:
                        val = Math.sin(c*Xi[2]);
                        break;
                    case 5:
                        val = 0;
                        break;
                    case 6:
                        val = 0;
                        break;
                    case 7:
                        val = Math.cos(c*Xi[2]);
                        break;
                    case 8:
                        val = 0;
                        break;
                }
                return val;
            case 1:
                switch(m) {
                    case 0:
                        val = 0;
                        break;
                    case 1:
                        val = 0;
                        break;
                    case 2:
                        val = 1;
                        break;
                    case 3:
                        val = 0;
                        break;
                    case 4:
                        val = 0;
                        break;
                    case 5:
                        val = Math.sin(c*Xi[2]);
                        break;
                    case 6:
                        val = 0;
                        break;
                    case 7:
                        val = 0;
                        break;
                    case 8:
                        val = Math.cos(c*Xi[2]);
                        break;
                }
                return val;
            case 2:
                
                switch(m) {
                    case 1:
                        val = 0;
                        break;
                    case 2:
                        val = 0;
                        break;
                    case 3:
                        val = c*Math.cos(c*Xi[2]);
                        break;
                    case 4:
                        val = c*Xi[0]*Math.cos(c*Xi[2]);
                        break;
                    case 5:
                        val = c*Xi[1]*Math.cos(c*Xi[2]);
                        break;
                    case 6:
                        val = -c*Math.sin(c*Xi[2]);
                        break;
                    case 7:
                        val = -c*Xi[0]*Math.sin(c*Xi[2]);
                        break;
                    case 8:
                        val = -c*Xi[1]*Math.sin(c*Xi[2]);
                        break;
                }
                return val;
        }
        return val;
    }
}

