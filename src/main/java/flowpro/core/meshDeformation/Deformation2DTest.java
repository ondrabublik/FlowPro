
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
public class Deformation2DTest extends Deformation {

    private final FluidForces[] fluidForces;
    private double t = 0;
    
    public Deformation2DTest(Parameters par, Equation eqn, int[][] TEale) throws IOException {
        super(par,eqn,TEale);
		
		fluidForces = new FluidForces[nBodies];
    }

	@Override
    public void newMeshPositionAndVelocity(Element[] elems, int timeOrder, double dt, double dto, MeshMove[] mshMov) {
        
        if(dt > 1000){
            dt = 0;
        }
        
        double a1 = 1.0 / dt;
        double a2 = -1.0 / dt;
        double a3 = 0;

        if (timeOrder > 1) {
            a1 = (2 * dt + dto) / (dt * (dt + dto));  // 3/(2*dt); 
            a2 = -(dt + dto) / (dt * dto);  // -2/dt;
            a3 = dt / (dto * (dt + dto));  // 1/(2*dt);
        }
        
        t = t + dt;
        // multiple blending function
        for (int i = 0; i < elems.length; i++) {
            for (int j = 0; j < elems[i].nVertices; j++) {
                
                //________________ test
                double L = 20;
                double t0 = 10; //40.1978;
                
                double A = 2;
                elems[i].vertices[j][0] = elems[i].vertices0[j][0] + A * Math.sin(2 * Math.PI * t/t0) * Math.sin(2 * Math.PI * elems[i].vertices0[j][0]/L) * Math.sin(2 * Math.PI * elems[i].vertices0[j][1]/L);
                elems[i].vertices[j][1] = elems[i].vertices0[j][1] + A * Math.sin(2 * Math.PI * t/t0) * Math.sin(2 * Math.PI * elems[i].vertices0[j][0]/L) * Math.sin(2 * Math.PI * elems[i].vertices0[j][1]/L);
                //________________ end test
                
                for (int d = 0; d < elems[i].dim; d++) {
                    elems[i].U[j][d] = a1 * elems[i].vertices[j][d] + a2 * elems[i].verticesOld[j][d] + a3 * elems[i].verticesOld2[j][d];
                }
            }
        }
    }

	@Override
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
    
	@Override
    public void calculateForces(Element[] elems, MeshMove[] mshMov) {        
        for (int b = 0; b < nBodies; b++) {
			double[] totalTranslationForce = new double[2];
			double[] totalRotationForce = new double[1];
            for (Element elem : elems) {
                for (int k = 0; k < elem.nFaces; k++) {
                    if (elem.TEale[k] == b + 2 && elem.insideMetisDomain) {
                        double[] Jac = elem.Int.faces[k].JacobianFace;
                        double[] weights = elem.Int.faces[k].weightsFace;
                        double[][] baseLeft = elem.Int.faces[k].basisFaceLeft;
                        double fx = 0;
                        double fy = 0;
                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
                            double[] WL = new double[elem.getNEqs()];
                            for (int j = 0; j < elem.getNEqs(); j++) {
                                for (int m = 0; m < elem.nBasis; m++) {
                                    WL[j] = WL[j] + elem.W[j * elem.nBasis + m] * baseLeft[p][m];
                                }
                            }
                            double pressure = eqn.pressure(WL);
                            fx += Jac[p] * weights[p] * elem.n[k][p][0] * pressure;
                            fy += Jac[p] * weights[p] * elem.n[k][p][1] * pressure;
                        }
                        totalTranslationForce[0] += fx;
                        totalTranslationForce[1] += fy;
                        totalRotationForce[0] += -fx * (elem.Xes[k][1] - center[1][b]) + fy * (elem.Xes[k][0] - center[0][b]);
                    }
                }
            }
			fluidForces[b] = new FluidForces(totalTranslationForce, totalRotationForce);
        }

        // user defined totalTranslationForce term
//        userDef = new double[nUserValues][nBodies];
//        for (int b = 0; b < nBodies; b++) {
//            for (Element elem : elems) {
//                for (int k = 0; k < elem.nFaces; k++) {
//                    if (elem.TEale[k] == b + 2 && elem.insideMetisDomain) {
//                        double[] Jac = elem.Int.faces[k].JacobianFace;
//                        double[] weights = elem.Int.faces[k].weightsFace;
//                        double[][] baseLeft = elem.Int.faces[k].basisFaceLeft;
//                        double fyTop = 0;
//                        double fyBottom = 0;
//                        for (int p = 0; p < elem.Int.faces[k].nIntEdge; p++) { // edge integral
//                            double[] WL = new double[elem.getNEqs()];
//                            for (int j = 0; j < elem.getNEqs(); j++) {
//                                for (int m = 0; m < elem.nBasis; m++) {
//                                    WL[j] = WL[j] + elem.W[j * elem.nBasis + m] * baseLeft[p][m];
//                                }
//                            }
//                            double pressure = eqn.pressure(WL);
//                            if (elem.n[k][p][1] > 0) {
//                                fyTop += Jac[p] * weights[p] * elem.n[k][p][1] * pressure;
//                            } else {
//                                fyBottom += Jac[p] * weights[p] * elem.n[k][p][1] * pressure;
//                            }
//                        }
//                        userDef[0][b] += fyTop;
//                        userDef[1][b] += fyBottom;
//                    }
//                }
//            }
//        }
    }
    
	@Override
    public FluidForces[] getFluidForces(){
        return fluidForces;
    }
}
