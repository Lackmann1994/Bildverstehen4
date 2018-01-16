import Jama.Matrix;
import Jama.QRDecomposition;
import Jama.SingularValueDecomposition;

import java.util.HashMap;
import java.util.Map;

public class Main {
    public static Point P1 = new Point(185.0, 49.0);
    public static Point P2 = new Point(503.0, 137.0);
    public static Point P3 = new Point(518.0, 434.0);
    public static Point P4 = new Point(163.0, 351.0);
    public static Point P1n = new Point(200.0, 200.0);
    public static Point P2n = new Point(600.0, 200.0);
    public static Point P3n = new Point(600.0, 600.0);
    public static Point P4n = new Point(200.0, 600.0);

    public static void main(String[] args) {
        System.out.println("Aufgabe 1a1:");
        Matrix H = solveInhomogenousLGS();
        H.print(10, 9);


        System.out.println("Aufgabe 1a2:");
        Matrix H2 = solveHomogenousLGS();
        H2.print(10, 9);


        System.out.println("Aufgabe 1b:");
        Map<String, Matrix> matrices = decomposeHomography(H);


        for (Map.Entry<String, Matrix> e : matrices.entrySet()) {
            System.out.println(e.getKey());
            e.getValue().print(10, 9);
        }

        //Aufgabe 2
        System.out.println("Aufgabe 2:");
        Matrix mP = new Matrix(new double[][]{
                {-926.70775, -925.63171, 1314.4708, 1423086.2},
                {1206.8268, -428.23672, 1106.7185, 524738.83},
                {-0.11055, 0.35751, 0.92734, 1001.7}});
        Matrix mnP = new Matrix(new double[][]{
                {-925.60853, 1208.7902, 1012.6846, 768151.47},
                {-58.627242, -366.6402, 1146.5259, 644372.57},
                {0.36806, 0.43941, 0.81941, 801.5}});
        Matrix F = fundamentalmatrix(mP, mnP);
        //Aufgabe 2.1
        System.out.println("Aufgabe2.2");
        System.out.println("fundamentalmatrix(mP, mnP)");
        F.print(10, 10);
        System.out.println("fundamentalmatrix(mnP, mP)");
        fundamentalmatrix(mnP, mP).transpose().print(10, 10);

        //Aufgabe 2.2
        System.out.println("Aufgabe2.2");
        Matrix point = new Matrix(new double[][]{
                {47},
                {11},
                {42},
                {1.0}
        });
        mnP.times(point).transpose().times(F).times(mP.times(point)).print(10, 10);
    }

    //Aufgabe 1a1
    public static Matrix solveInhomogenousLGS() {
        Matrix A = new Matrix(new double[][]
                {{P1.X, P1.Y, 0.0, 0.0, 1.0, 0.0, P1.X * P1n.X, P1.Y * P1n.X},
                        {P2.X, P2.Y, 0.0, 0.0, 1.0, 0.0, P2.X * P2n.X, P2.Y * P2n.X},
                        {P3.X, P3.Y, 0.0, 0.0, 1.0, 0.0, P3.X * P3n.X, P3.Y * P3n.X},
                        {P4.X, P4.Y, 0.0, 0.0, 1.0, 0.0, P4.X * P4n.X, P4.Y * P4n.X},
                        {0.0, 0.0, P1.X, P1.Y, 0.0, 1.0, P1.X * P1n.Y, P1.Y * P1n.Y},
                        {0.0, 0.0, P2.X, P2.Y, 0.0, 1.0, P2.X * P2n.Y, P2.Y * P2n.Y},
                        {0.0, 0.0, P3.X, P3.Y, 0.0, 1.0, P3.X * P3n.Y, P3.Y * P3n.Y},
                        {0.0, 0.0, P4.X, P4.Y, 0.0, 1.0, P4.X * P4n.Y, P4.Y * P4n.Y}
                });
        Matrix b = new Matrix(new double[][]{{P1n.X}, {P2n.X}, {P3n.X}, {P4n.X}, {P1n.Y}, {P2n.Y}, {P3n.Y}, {P4n.Y}});
        Matrix x = A.solve(b);


        double[][] matH = {{x.get(0, 0), x.get(1, 0), x.get(4, 0)},
                {x.get(2, 0), x.get(3, 0), x.get(5, 0)},
                {x.get(6, 0), x.get(7, 0), 1.0}};
        return new Matrix(matH);
    }

    //Aufgabe 1a2
    public static Matrix solveHomogenousLGS() {
        Matrix A = new Matrix(new double[][]
                {{0.0, 0.0, 0.0, -P1.X, -P1.Y, 1.0, P1.X * P1n.Y, P1.Y * P1n.Y, P1n.Y},
                        {0.0, 0.0, 0.0, -P2.X, -P2.Y, 1.0, P2.X * P2n.Y, P2.Y * P2n.Y, P2n.Y},
                        {0.0, 0.0, 0.0, -P3.X, -P3.Y, 1.0, P3.X * P3n.Y, P3.Y * P3n.Y, P3n.Y},
                        {0.0, 0.0, 0.0, -P4.X, -P4.Y, 1.0, P4.X * P4n.Y, P4.Y * P4n.Y, P4n.Y},
                        {P1.X, P1.Y, 1.0, 0.0, 0.0, 0.0, -P1.X * P1n.X, -P1.Y * P1n.X, -P1n.X},
                        {P2.X, P2.Y, 1.0, 0.0, 0.0, 0.0, -P2.X * P2n.X, -P2.Y * P2n.X, -P2n.X},
                        {P3.X, P3.Y, 1.0, 0.0, 0.0, 0.0, -P3.X * P3n.X, -P3.Y * P3n.X, -P3n.X},
                        {P4.X, P4.Y, 1.0, 0.0, 0.0, 0.0, -P4.X * P4n.X, -P4.Y * P4n.X, -P4n.X},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                });

        Matrix V = new SingularValueDecomposition(A).getV();
        double factor = 1 / V.get(8, 8);
        return new Matrix(new double[][]{{V.get(0, 8) * factor, V.get(1, 8) * factor, V.get(2, 8) * factor},
                {V.get(3, 8) * factor, V.get(4, 8) * factor, V.get(5, 8) * factor},
                {V.get(6, 8) * factor, V.get(7, 8) * factor, 1.0}});
    }


    //Aufgabe 1b
    public static Map<String, Matrix> decomposeHomography(Matrix H) {

        Matrix A = H.getMatrix(0, 1, 0, 1);
        Matrix tv = H.getMatrix(0, 1, 2, 2).times(H.getMatrix(2, 2, 0, 1));
        double p3 = H.get(2, 2);

        QRDecomposition QRD = new QRDecomposition(A.minus(tv));

        Matrix sR = QRD.getQ();
        Matrix K = QRD.getR();

        double factor = Math.sqrt(1 / (K.get(0, 0) * K.get(1, 1)));

        sR.timesEquals(1 / factor);
        K.timesEquals(factor);

        Matrix Hs = new Matrix(3, 3);
        Hs.setMatrix(0, 1, 0, 1, sR);
        Hs.setMatrix(0, 1, 2, 2, H.getMatrix(0, 1, 2, 2));
        Hs.set(2, 2, 1.0);

        Matrix Ha = new Matrix(3, 3);
        Ha.setMatrix(0, 1, 0, 1, K);
        Ha.set(2, 2, 1.0);

        Matrix Hp = new Matrix(3, 3);
        Hp.set(0, 0, 1.0);
        Hp.set(1, 1, 1.0);
        Hp.set(2, 2, p3);
        Hp.setMatrix(2, 2, 0, 1, H.getMatrix(2, 2, 0, 1));

        Map<String, Matrix> matrices = new HashMap<>();
        matrices.put("sR", sR);
        matrices.put("K", K);
        matrices.put("H", H);
        matrices.put("Ha", Ha);
        matrices.put("Hp", Hp);
        matrices.put("Hs", Hs);
        return matrices;
    }


    //Aufgabe 2
    public static Matrix fundamentalmatrix(Matrix P, Matrix Pn) {
        // Pseudo Rechtsinverse
        Matrix rightPseudoInverse = P.transpose().times(P.times(P.transpose()).inverse());


        //Berechnung von C (siehe Skript 37)
        Matrix X = P.getMatrix(0, 2, 1, 3);
        Matrix Y = new Matrix(3, 3);
        Y.setMatrix(0, 2, 0, 0, P.getMatrix(0, 2, 0, 0));
        Y.setMatrix(0, 2, 1, 2, P.getMatrix(0, 2, 2, 3));
        Matrix Z = new Matrix(3, 3);
        Z.setMatrix(0, 2, 0, 1, P.getMatrix(0, 2, 0, 1));
        Z.setMatrix(0, 2, 2, 2, P.getMatrix(0, 2, 3, 3));
        Matrix H = P.getMatrix(0, 2, 0, 2);
        Matrix C = new Matrix(new double[][]{
                {X.det() / -H.det()},
                {-Y.det() / -H.det()},
                {Z.det() / -H.det()},
                {1.0}});


        //Berechnung von [e´]x (Siehe Skript 86)
        Matrix epipole = Pn.times(C);
        epipole.set(0, 0, epipole.get(0, 0) / epipole.get(2, 0));
        epipole.set(1, 0, epipole.get(1, 0) / epipole.get(2, 0));
        epipole.set(2, 0, 1);

        Matrix epipoleX = new Matrix(new double[][]{
                {0.0, -epipole.get(2, 0), epipole.get(1, 0)},
                {epipole.get(2, 0), 0.0, -epipole.get(0, 0)},
                {-epipole.get(1, 0), epipole.get(0, 0), 0.0}});


        Matrix foundamental = epipoleX.times(Pn.times(rightPseudoInverse));

        double factor = foundamental.get(2, 2);

        for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
                foundamental.set(ii, jj, foundamental.get(ii, jj) / factor);
            }
        }


//        Map<String, Matrix> matrices = new HashMap<>();
//        matrices.put("C", C);
//        matrices.put("epipole", epipole);
//        matrices.put("foundamental", foundamental);
//        for (Map.Entry<String, Matrix> e : matrices.entrySet()) {
//            System.out.println(e.getKey());
//            e.getValue().print(10, 9);
//        }

        return foundamental;
    }
}


class Point {
    public double X;
    public double Y;
    public Point(double x, double y) {
        X = x;
        Y = y;
    }
}
