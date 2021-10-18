package flowpro.core.curvedBoundary;

import flowpro.api.Mat;
import flowpro.core.elementType.ElementType;
import java.io.IOException;
import static flowpro.core.elementType.ElementType.firstDigit;

/**
 *
 * @author obublik
 */
public class CurvedBoundary {

    public static FaceCurvature[] modifyMesh(int[] elemsType, double[][] PXY, int[][] TP, int[][] TT) {

        FaceCurvature[] fCurv = new FaceCurvature[TP.length];
        elemsType = ElementType.firstDigit(elemsType); // reset curved elements
        
        // compute point normals
        double[][] t = new double[PXY.length][];
        for (int i = 0; i < TP.length; i++) {
            if (ElementType.firstDigit(elemsType[i]) == 3) {
                for (int j = 0; j < TT[i].length; j++) {
                    if (TT[i][j] == -1) {
                        double[] tangent = new double[]{PXY[TP[i][(j + 1) % 3]][0] - PXY[TP[i][j]][0], PXY[TP[i][(j + 1) % 3]][1] - PXY[TP[i][j]][1]};
                        tangent = Mat.normVector(tangent);
                        if (t[TP[i][j]] == null) {
                            t[TP[i][j]] = new double[]{tangent[0], tangent[1]};
                        } else {
                            t[TP[i][j]][0] = (t[TP[i][j]][0] + tangent[0]) / 2;
                            t[TP[i][j]][1] = (t[TP[i][j]][1] + tangent[1]) / 2;
                        }
                        int jp = (j + 1) % 3;
                        if (t[TP[i][jp]] == null) {
                            t[TP[i][jp]] = new double[]{tangent[0], tangent[1]};
                        } else {
                            t[TP[i][jp]][0] = (t[TP[i][jp]][0] + tangent[0]) / 2;
                            t[TP[i][jp]][1] = (t[TP[i][jp]][1] + tangent[1]) / 2;
                        }
                        break;
                    }
                }
            }
        }

        // find curved boundary
        double tolTop = 1 - 1e-12;
        double tolBottom = 0.1;
        for (int i = 0; i < TP.length; i++) {
            if (firstDigit(elemsType[i]) == 3) {
                for (int j = 0; j < TT[i].length; j++) {
                    if (TT[i][j] == -1) {
                        double[] tangent = new double[]{PXY[TP[i][(j + 1) % 3]][0] - PXY[TP[i][j]][0], PXY[TP[i][(j + 1) % 3]][1] - PXY[TP[i][j]][1]};
                        tangent = Mat.normVector(tangent);
                        double scal1 = Math.abs(Mat.scalar(tangent, t[TP[i][j]]));
                        double scal2 = Math.abs(Mat.scalar(tangent, t[TP[i][(j + 1) % 3]]));
                        if ((scal1 < tolTop && scal1 > tolBottom) && (scal2 < tolTop && scal2 > tolBottom)) {
                            fCurv[i] = new Curved2DLine(t[TP[i][j]], t[TP[i][(j + 1) % 3]], PXY[TP[i][j]], PXY[TP[i][(j + 1) % 3]]);
                            // permute vertices, change element type
                            permuteTriangle(j, i, TT, TP);
                            elemsType[i] = 31;
                        } else {
                            if ((scal1 >= tolTop && scal2 < tolTop) || (scal1 <= tolBottom && scal2 > tolBottom)) {
                                fCurv[i] = new Curved2DLine("right",t[TP[i][(j + 1) % 3]], PXY[TP[i][j]], PXY[TP[i][(j + 1) % 3]]);
                                // permute vertices, change element type
                                permuteTriangle(j, i, TT, TP);
                                elemsType[i] = 31;
                            } else {
                                if((scal1 < tolTop && scal2 >= tolTop) || (scal1 > tolBottom && scal2 <= tolBottom)){
                                    fCurv[i] = new Curved2DLine("left",t[TP[i][j]], PXY[TP[i][j]], PXY[TP[i][(j + 1) % 3]]);
                                // permute vertices, change element type
                                permuteTriangle(j, i, TT, TP);
                                elemsType[i] = 31;
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }

        return fCurv;
    }

    public static void saveMesh(String geometryPath, int[] elemsType, int[][] TP, int[][] TT) throws IOException {
        Mat.save(elemsType, geometryPath + "typ.txt");
        Mat.save(TP, geometryPath + "TP.txt");
        Mat.save(TT, geometryPath + "TT.txt");
    }

    private static void permuteTriangle(int j, int i, int[][] TT, int[][] TP) {
        if (j == 0) {
            TP[i] = new int[]{TP[i][1], TP[i][2], TP[i][0]};
            TT[i] = new int[]{TT[i][1], TT[i][2], TT[i][0]};
        }
        if (j == 1) {
            TP[i] = new int[]{TP[i][2], TP[i][0], TP[i][1]};
            TT[i] = new int[]{TT[i][2], TT[i][0], TT[i][1]};
        }
    }
}
