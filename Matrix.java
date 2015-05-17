/*************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones collection of static methods for manipulating
 *  matrices.
 *
 *************************************************************************/

import java.math.*;
import org.apache.commons.math3.fraction.BigFraction;

public class Matrix {
    
    public static BigFraction[][] clone(BigFraction[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length;
        BigFraction[][] clone = new BigFraction[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                clone[i][j] = new BigFraction(matrix[i][j].getNumerator(), matrix[i][j].getDenominator());
        return clone;
    }
    
    // return a random m-by-n matrix with values between 0 and 1
    public static BigFraction[][] random(int m, int n) {
        BigFraction[][] C = new BigFraction[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = new BigFraction(Math.random());
        return C;
    }
    
    // return a random m-by-n matrix with values between 0 and 10
    public static BigFraction[][] randomInt(int m, int n) {
        BigFraction[][] C = new BigFraction[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = new BigFraction(Math.round(Math.random() * 20) - 10);// * 256;
        return C;
    }
    
    // return n-by-n identity matrix I
    public static BigFraction[][] identity(int n) {
        BigFraction[][] I = new BigFraction[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                I[i][j] = new BigFraction(0);
                I[i][i] = new BigFraction(1);
            }
        }
        return I;
    }
    
    // return x^T y
    public static BigFraction dot(BigFraction[] x, BigFraction[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        BigFraction sum = new BigFraction(0);
        for (int i = 0; i < x.length; i++)
            sum.add(x[i].multiply(y[i]));
        return sum;
    }
    
    // return C = A^T
    public static BigFraction[][] transpose(BigFraction[][] A) {
        int m = A.length;
        int n = A[0].length;
        BigFraction[][] C = new BigFraction[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[j][i] = A[i][j];
        return C;
    }
    
    // return C = A + B
    public static BigFraction[][] add(BigFraction[][] A, BigFraction[][] B) {
        int m = A.length;
        int n = A[0].length;
        BigFraction[][] C = new BigFraction[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j].add(B[i][j]);
        return C;
    }
    
    // return C = A - B
    public static BigFraction[][] subtract(BigFraction[][] A, BigFraction[][] B) {
        int m = A.length;
        int n = A[0].length;
        BigFraction[][] C = new BigFraction[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j].subtract(B[i][j]);
        return C;
    }
    
    // return C = A * B
    public static BigFraction[][] multiply(BigFraction[][] A, BigFraction[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (nA != mB) throw new RuntimeException("Illegal matrix dimensions.");
        BigFraction[][] C = new BigFraction[mA][nB];
        for (int i = 0; i < mA; i++)
            for (int j = 0; j < nB; j++)
                for (int k = 0; k < nA; k++)
                    C[i][j] = (A[i][k].multiply(B[k][j]));
        return C;
    }
    
    // matrix-vector multiplication (y = A * x)
    public static BigFraction[] multiply(BigFraction[][] A, BigFraction[] x) {
        int m = A.length;
        int n = A[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        BigFraction[] y = new BigFraction[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i].add(A[i][j].multiply(x[j]));
        return y;
    }
    
    
    // vector-matrix multiplication (y = x^T A)
    public static BigFraction[] multiply(BigFraction[] x, BigFraction[][] A) {
        int m = A.length;
        int n = A[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        BigFraction[] y = new BigFraction[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j].add(A[i][j].multiply(x[i]));
        return y;
    }
    
    
    // matrix-scalar multiplication (B = Ac)
    public static BigFraction[][] multiply(BigFraction[][] A, BigFraction c) {
        int m = A.length;
        int n = A[0].length;
        BigFraction[][] B = new BigFraction[n][m];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                B[i][j] = c.multiply(A[i][j]);
        return B;
    }
    
    
    //vector-scalar multiplication (y = vector * scalar)
    public static BigFraction[] vectorTimesScalar(BigFraction[] vector, BigFraction scalar) {
        vector = vector.clone();
        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i].multiply(scalar);
        }
        return vector;
    }
    
    
    //vector-vector addition (y = a + b)
    public static BigFraction[] vectorAddition(BigFraction[] a, BigFraction[] b) {
        a = a.clone();
        b = b.clone();
        if (a.length != b.length) throw new RuntimeException("Illegal vector dimensions.");
        for (int i = 0; i < a.length; i++) {
            a[i] = a[i].add(b[i]);
        }
        return a;
    }
    
    // test client
    public static void main(String[] args) {
        System.out.println("D");
        System.out.println("--------------------");
        BigFraction[][] d = { { new BigFraction(1), new BigFraction(2), new BigFraction(3) },
            { new BigFraction(4), new BigFraction(5), new BigFraction(6) },
            { new BigFraction(7), new BigFraction(8), new BigFraction(9)} };
        LinearAlgebra.print(d);
        System.out.println();
        
        System.out.println("I");
        System.out.println("--------------------");
        BigFraction[][] c = Matrix.identity(5);
        LinearAlgebra.print(c);
        System.out.println();
        
        System.out.println("A");
        System.out.println("--------------------");
        BigFraction[][] a = Matrix.random(5, 5);
        LinearAlgebra.print(a);
        System.out.println();
        
        System.out.println("A^T");
        System.out.println("--------------------");
        BigFraction[][] b = Matrix.transpose(a);
        LinearAlgebra.print(b);
        System.out.println();
        
        System.out.println("A + A^T");
        System.out.println("--------------------");
        BigFraction[][] e = Matrix.add(a, b);
        LinearAlgebra.print(e);
        System.out.println();
        
        System.out.println("A * A^T");
        System.out.println("--------------------");
        BigFraction[][] f = Matrix.multiply(a, b);
        LinearAlgebra.print(f);
        System.out.println();
        
        System.out.println("c * A");
        System.out.println("--------------------");
        BigFraction scalar = new BigFraction(2);
        BigFraction[][] g = Matrix.multiply(a, scalar);
        LinearAlgebra.print(g);
        System.out.println();
        
        System.out.println("c * v");
        System.out.println("--------------------");
        BigFraction[] v = { new BigFraction(1), new BigFraction(2), new BigFraction(3) };
        BigFraction scalar2 = new BigFraction(3);
        BigFraction[] h = Matrix.vectorTimesScalar(v, scalar2);
        LinearAlgebra.print(h);
        System.out.println();
        
        System.out.println("v + u");
        System.out.println("--------------------");
        BigFraction[] v2 = { new BigFraction(1), new BigFraction(2), new BigFraction(3) };
        BigFraction[] u = { new BigFraction(9), new BigFraction(8), new BigFraction(7) };
        BigFraction[] i = Matrix.vectorAddition(v2, u);
        LinearAlgebra.print(i);
        System.out.println();
    }
}