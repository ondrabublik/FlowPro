package flowpro.core.meshDeformation;

import flowpro.api.Equation;
import flowpro.api.FluidForces;
import flowpro.api.MeshMove;
import flowpro.core.element.Element;
import flowpro.core.Parameters;
import java.io.*;

/**
 *
 * @author obublik
 */
public class Deformation3DBeam extends Deformation {

    public double[][] boundaryForce;
    public int[][] faceIndexes;

    public Deformation3DBeam(Parameters par, Equation eqn, int[][] TEale) throws IOException {
        super(par, eqn, TEale);
    }

    public void newMeshPositionAndVelocity(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {

        double a1 = 1.0 / dt;
        double a2 = -1.0 / dt;
        double a3 = 0;

        if (timeOrder > 1) {
            a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
            a2 = -(dt + dto) / (dt * dto);  // -2/dt;
            a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
        }

        // multiple blending function
        for (int i = 0; i < elems.length; i++) {
            for (int j = 0; j < elems[i].nVertices; j++) {
                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].vertices[j][d] = 0;
                }
                
                for (int k = 0; k < nBodies; k++) {
                    double[][] rbfCoeffs = mshMov[k].getDeformationCoefficients();
                    double xTrans = 0;
                    double yTrans = 0;
                    double alpha = 0;
                    if (rbfCoeffs != null) {
                        double[][] beamZcoord = mshMov[k].getBoundaryPoints();
                        for (int p = 0; p < rbfCoeffs[0].length; p++) {
                            double rbf = radialBasisFunction(elems[i].vertices0[j][2], beamZcoord[p][0]);
                            xTrans += rbfCoeffs[0][p] * rbf;
                            yTrans += rbfCoeffs[1][p] * rbf;
                            alpha += rbfCoeffs[2][p] * rbf;
                        }
                    }
                    
                    double xNew = (Math.cos(alpha) * (elems[i].vertices0[j][0] - center[0][k]) - Math.sin(alpha) * (elems[i].vertices0[j][1] - center[0][k]) + center[0][k] - elems[i].vertices0[j][0]) * nBodies + elems[i].vertices0[j][0] + xTrans * nBodies;
                    double yNew = (Math.sin(alpha) * (elems[i].vertices0[j][0] - center[0][k]) + Math.cos(alpha) * (elems[i].vertices0[j][1] - center[1][k]) + center[1][k] - elems[i].vertices0[j][1]) * nBodies + elems[i].vertices0[j][1] + yTrans * nBodies;               
                    elems[i].vertices[j][0] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][0] + elems[i].blendFun[j][k] * xNew) / nBodies;
                    elems[i].vertices[j][1] += ((1 - elems[i].blendFun[j][k]) * elems[i].vertices0[j][1] + elems[i].blendFun[j][k] * yNew) / nBodies;
                }
                elems[i].vertices[j][2] = elems[i].vertices0[j][2];

                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].U[j][d] = a1 * elems[i].vertices[j][d] + a2 * elems[i].verticesOld[j][d] + a3 * elems[i].verticesOld2[j][d];
                }
            }
        }
    }

    public void nextTimeLevel(Element[] elems) {
        for (int i = 0; i < elems.length; i++) {
            for (int j = 0; j < elems[i].nVertices; j++) {
                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].verticesOld2[j][d] = elems[i].verticesOld[j][d];
                    elems[i].verticesOld[j][d] = elems[i].vertices[j][d];
                }
            }
        }
    }

    public void calculateForces(Element[] elems, MeshMove[] mshMov) {
        boundaryForce = null;
        faceIndexes = null;

        int s = 0;
        for (Element elem : elems) {
            for (int k = 0; k < elem.nFaces; k++) {
                if (elem.TEale[k] > 1 && elem.insideMetisDomain) {
                    s++;
                }

            }
        }
        
        if (s > 0) {
            int dim = elems[0].dim;
            boundaryForce = new double[s][2*dim]; // s x [Fx,Fy,Fz,xes,yes,zes]
            faceIndexes = new int[s][3];
            s = 0;
            for (int i = 0; i < elems.length; i++) {
                Element elem = elems[i];
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TEale[k] > 1 && elem.insideMetisDomain) {
                        double pressure = eqn.pressure(elem.calculateAvgW());
                        for (int d = 0; d < dim; d++) {
                            boundaryForce[s][d] = pressure;
                        }
                        System.arraycopy(elem.Xes[k], 0, boundaryForce[s], dim, dim);
                        faceIndexes[s][0] = i;
                        faceIndexes[s][1] = k;
                        faceIndexes[s][2] = elem.TEale[k];
                        s++;
                    }
                }
            }
        }
    }

    public FluidForces[] getFluidForces() {
        return null;  // prozatim abych prelozil kod
    }

    public double radialBasisFunction(double a, double b) {
        return Math.max(1 - Math.abs(a - b),0);
    }
}
