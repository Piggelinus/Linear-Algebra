import java.util.*;
import java.math.*;
import org.apache.commons.math3.fraction.BigFraction;

public class LinearAlgebra extends Matrix {
    
    private BigFraction[][] matrix;
    private List<BigFraction[][]> savedMatrices = new ArrayList<BigFraction[][]>();
    
    public static void main(String[] args) {
        try {
            LinearAlgebra la;
            if (args.length == 2) {
                la = new LinearAlgebra(new Integer(args[0]), new Integer(args[1]));
            }
            else {
                la = new LinearAlgebra(3, 3);
            }
            la.manage();
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }
    
    public LinearAlgebra(int m, int n) {
        matrix = randomInt(m, n);
        savedMatrices.add(Matrix.clone(matrix));
    }
    
    public BigFraction[][] getMatrix() {
        return matrix;
    }
    
    public BigFraction[][] newCustomMatrix(int n) {
        Scanner s = new Scanner(System.in);
        BigFraction[][] A = new BigFraction[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = new BigFraction(s.nextInt());
            }
        }
        return A;
    }
    
    public BigFraction[][] newCustomMatrix(int n, int m) {
        Scanner s = new Scanner(System.in);
        BigFraction[][] A = new BigFraction[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = new BigFraction(s.nextInt());
            }
        }
        return A;
    }
    
    public BigFraction[] customVector(int n) {
        Scanner s = new Scanner(System.in);
        BigFraction[] X = new BigFraction[n];
        for (int i = 0; i < n; i++) {
            X[i] = new BigFraction(s.nextInt());
        }
        return X;
    }
    
    public static void print(BigFraction[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                System.out.print(matrix[i][j].toString().replaceAll(" ", "") + " ");
            }
            System.out.println();
        }
    }
    
    public static void print(BigFraction[] matrix) {
        int n = matrix.length;
        for (int i = 0; i < n; i++) {
            System.out.print(matrix[i].toString().replaceAll(" ", "") + " ");
        }
        System.out.println();
    }
    
    public void help() {
        System.out.println("Help - Listing all available commands");
        System.out.println("-------------------------------------");
        System.out.println("Matrices and Vectors:");
        System.out.println("Elementary row operations:");
        System.out.println("swap a b : swaps row a and b (index starts at 0)");
        System.out.println("multiply/*/x row scalar : multiplies row a with scalar");
        System.out.println("add/+ a b scalar : adds a multiple of row a (a * scalar) to row b");
        System.out.println("set : multiplies a row by a factor to set the leading variable to 1");
        System.out.println("det/d : calculates the determinant");
        System.out.println("trace/tr : calculates the trace");
        System.out.println("transpose/t : transposes the matrix");
        System.out.println("inverse/i : inverts the matrix");
        System.out.println("ref : coverts matrix into reduced row echelon form");
        System.out.println("solve : finds the basis of the null space");
        System.out.println("eigenvalues/ev : finds the eigenvalues of the matrix");
        System.out.println("eigenvectorbasis/evb : finds eigenvalues and the eigenvectors that form the basis for each eigenspace");
        System.out.println("-------------------------------------");
        System.out.println("Polynomials:");
        System.out.println("q a b c : find roots for quadratic, given coefficients a, b & c");
        System.out.println("-------------------------------------");
    }
    
    public void manage() {
        System.out.println("Here is your matrix");
        print(getMatrix());
        Scanner scanner = new Scanner(System.in);
        while (scanner.hasNext()) {
            String input = scanner.next();
            if (input.toLowerCase().equals("help") || input.toLowerCase().equals("h")) {
                help();
            }
            if (input.toLowerCase().equals("swap")) {
                try {
                    swapRows(matrix, new Integer(scanner.next()) - 1, new Integer(scanner.next()) - 1);
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("multiply") || input.equals("*") || input.equals("x")) {
                try {
                    multiplyRow(matrix, scanner.nextInt() - 1, new BigFraction(scanner.nextInt()));
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("multiply*") || input.equals("**") || input.equals("x*")) {
                try {
                    int scalar = scanner.nextInt();
                    for (int i = 0; i < matrix.length; i++) {
                        multiplyRow(matrix, i, new BigFraction(scalar));
                    }
                    
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("divide") || input.equals("/")) {
                try {
                    multiplyRow(matrix, scanner.nextInt() - 1, new BigFraction(1, scanner.nextInt()));
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("divide*") || input.equals("/*")) {
                try {
                    int scalar = scanner.nextInt();
                    for (int i = 0; i < matrix.length; i++) {
                        multiplyRow(matrix, i, new BigFraction(1, scalar));
                    }
                    
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("add") || input.equals("+")) {
                try {
                    addRow(matrix, new Integer(scanner.next()) - 1, new Integer(scanner.next()) - 1, new BigFraction(scanner.nextInt()));
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("set") || input.equals("+")) {
                try {
                    leadingOne(matrix, new Integer(scanner.next()) - 1);
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("set*") || input.equals("+")) {
                try {
                    for (int i = 0; i < matrix.length; i++) {
                        leadingOne(matrix, i);
                    }
                } catch (Exception e) {
                    
                }
            }
            else if (input.toLowerCase().equals("ref")) {
                try {
                    ref(matrix, 0, 0);
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("solve")) {
                try {
                    matrix = solve(matrix);
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("det") || input.toLowerCase().equals("d")) {
                try {
                    System.out.println(determinant(matrix));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("trace") || input.toLowerCase().equals("tr")) {
                try {
                    System.out.println(trace(matrix));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            /*
            else if (input.toLowerCase().equals("augment") || input.toLowerCase().equals("au")) {
                int n = matrix.length;
                int m = matrix[0].length;
                try {
                    BigFraction[] column = new BigFraction[m];
                    for (int i = 0; i < column.length; i++) {
                        column[i] = new BigFraction(scanner.nextInt());
                    }
                    BigFraction[][] newMatrix = new BigFraction[n][m + 1];
                    newMatrix = augment(matrix, column);
                    print(newMatrix);
                    matrix = new BigFraction[n][m + 1];
                    matrix = newMatrix;
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }*/
            else if (input.toLowerCase().equals("transpose") || input.toLowerCase().equals("t")) {
                try {
                    matrix = transpose(matrix);
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("minors") || input.toLowerCase().equals("m")) {
                try {
                    matrix = matrixOfMinors(matrix);
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("cofactors") || input.toLowerCase().equals("cof")) {
                try {
                    matrix = matrixOfCofactors(matrix);
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("inverse") || input.toLowerCase().equals("i")) {
                try {
                    matrix = inverse(matrix);
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("load") || input.toLowerCase().equals("l")) {
                try {
                    matrix = Matrix.clone(savedMatrices.get(scanner.nextInt()));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("save") || input.toLowerCase().equals("s")) {
                try {
                    savedMatrices.set(scanner.nextInt(), Matrix.clone(matrix));
                } catch (Exception e) {
                    try {
                        savedMatrices.add(Matrix.clone(matrix));
                    } catch (Exception f) {
                        System.out.println("I don't know why, but you can't save non-square matrices :P");
                    }
                }
            }
            else if (input.toLowerCase().equals("remove") || input.toLowerCase().equals("r")) {
                try {
                    savedMatrices.remove(scanner.nextInt());
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("clear")) {
                try {
                    savedMatrices.clear();
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("list") || input.toLowerCase().equals("li")) {
                try {
                    System.out.println("-----------------------");
                    System.out.println("listing saved matrices:");
                    for (int i = 0; i < savedMatrices.size(); i++) {
                        System.out.println("index : " + i);
                        print(savedMatrices.get(i));
                        System.out.println();
                    }
                    System.out.println("-----------------------");
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("new")) {
                try {
                    matrix = randomInt(new Integer(scanner.next()), new Integer(scanner.next()));
                    savedMatrices.add(Matrix.clone(matrix));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("n")) {
                try {
                    matrix = newCustomMatrix(new Integer(scanner.next()));
                    savedMatrices.add(Matrix.clone(matrix));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("nm")) {
                try {
                    matrix = newCustomMatrix(new Integer(scanner.next()), new Integer(scanner.next()));
                    savedMatrices.add(Matrix.clone(matrix));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("eigenvalues") || input.toLowerCase().equals("ev")) {
                try {
                    System.out.print("Eigenvalues : ");
                    print(eigenValues(matrix));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("eigenvectorbasis") || input.toLowerCase().equals("evb")) {
                try {
                    BigFraction[] eigenValues = eigenValues(matrix);
                    for (int i = 0; i < eigenValues.length; i++) {
                        System.out.println("Basis for eigenspace corresponding to lambda = " + eigenValues[i]);
                        print(eigenVectorBasis(matrix.clone(), eigenValues[i]));
                    }
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("q")) {
                try {
                    print(quadraticSolver(new BigFraction(scanner.nextInt()), new BigFraction(scanner.nextInt()), new BigFraction(scanner.nextInt())));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("c")) {
                try {
                    System.out.println(combinations(new BigInteger(scanner.next()), new BigInteger(scanner.next())));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("p")) {
                try {
                    System.out.println(permutations(new BigInteger(scanner.next()), new BigInteger(scanner.next())));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("expand")) {
                try {
                    print(expandBinomial(new BigFraction(scanner.nextInt()), new BigFraction(scanner.nextInt()), new BigInteger(scanner.next())));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("expectedvalue")) {
                try {
                    BigFraction[] randomVariable = customVector(scanner.nextInt());
                    System.out.println(expectedValue(randomVariable));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            else if (input.toLowerCase().equals("variance")) {
                try {
                    BigFraction[] randomVariable = customVector(scanner.nextInt());
                    System.out.println(variance(randomVariable));
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
            print(getMatrix());
        }
    }
    
    public void swapRows(BigFraction[][] matrix, int a, int b) {
        BigFraction[] c = matrix[a];
        matrix[a] = matrix[b];
        matrix[b] = c;
    }
    
    public void multiplyRow(BigFraction[][] matrix, int row, BigFraction multiple) {
        matrix[row] = vectorTimesScalar(matrix[row], multiple);
    }
    
    public void addRow(BigFraction[][] matrix, int row, int destination, BigFraction multiple) {
        matrix[destination] = vectorAddition(matrix[destination], vectorTimesScalar(matrix[row], multiple));
    }
    
    public void leadingOne(BigFraction[][] matrix, int row) {
        BigFraction a = new BigFraction(0);
        for (int i = 0; i < matrix[row].length; i++) {
            if (!matrix[row][i].equals(BigFraction.ZERO)) {
                a = matrix[row][i];
                break;
            }
        }
        if (!a.equals(BigFraction.ZERO)) {
            BigFraction inverseOfa = BigFraction.ONE.divide(a);
            multiplyRow(matrix, row, inverseOfa);
        }
    }
    
    public void ref(BigFraction[][] matrix, int row, int column) {
        if (row < matrix.length && column < matrix.length) {
            for (int i = 0; i < matrix.length; i++) {
                leadingOne(matrix, i);
            }
            boolean allLeadingZeros = true;
            for (int i = row; i < matrix.length; i++) {
                if (matrix[i][column].equals(BigFraction.ONE)) {
                    allLeadingZeros = false;
                    if (i == row) {
                        break;
                    }
                    else {
                        swapRows(matrix, column, i);
                    }
                }
            }
            for (int i = row + 1; i < matrix.length; i++) {
                if (matrix[i][column].equals(BigFraction.ONE)) {
                    addRow(matrix, column, i, new BigFraction(-1));
                }
            }
            if (allLeadingZeros) {
                ref(matrix, row, column + 1);
            }
            else {
                ref(matrix, row + 1, column + 1);
            }
            for (int r = row - 1; r >= 0; r--) {
                addRow(matrix, row, r, matrix[r][column].negate());
            }
        }
    }
    
    public BigFraction[][] solve(BigFraction[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length;
        BigFraction[][] basis;
        ref(matrix, 0, 0);
        if (n == m && !determinant(matrix).equals(BigFraction.ZERO)) {
            BigFraction[] vector = new BigFraction[m];
            for (int i = 0; i < vector.length; i++) {
                vector[i] = BigFraction.ZERO;
            }
            basis = new BigFraction[1][m];
            basis[0] = vector;
        }
        else {
            Integer[] freeVariables = findFreeVariables(matrix);
            basis = new BigFraction[freeVariables.length][m];
            int num = 0;
            boolean looped = false;
            for (int i = 0; i < freeVariables.length; i++) {
                looped = true;
                BigFraction[] vector = new BigFraction[m];
                for (int row = 0; row < n; row++) {
                    vector[row] = BigFraction.ZERO;
                    for (int col = 0; col < m; col++) {
                        if (col == freeVariables[i]) {
                            vector[row] = vector[row].subtract(matrix[row][col]);
                        }
                    }
                    if (row == freeVariables[i]) {
                        vector[row] = BigFraction.ONE;
                    }
                }
                basis[i] = vector;
            }
        }
        return basis;
    }
    
    public BigFraction[][] augment(BigFraction[][] matrix, BigFraction[] column) {
        int n = matrix.length;
        int m = matrix[0].length;
        if (column.length != m)
            throw new RuntimeException("column to augment must have " + m + " columns");
        BigFraction[][] augmentedMatrix = new BigFraction[n][m + 1];
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < m + 1; col++) {
                if (col <= m - 1) {
                    augmentedMatrix[row][col] = matrix[row][col];
                }
                else {
                    augmentedMatrix[row][col] = column[col];
                }
            }
        }
        print(augmentedMatrix);
        return augmentedMatrix;
    }
    
    public Integer[] findFreeVariables(BigFraction[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length;
        ArrayList<Integer> freeVariables = new ArrayList<Integer>();
        for (int col = 0; col < m; col++) {
            ArrayList<BigFraction> values = new ArrayList<BigFraction>();
            for (int row = 0; row < n; row++) {
                if (matrix[row][col].compareTo(BigFraction.ZERO) != 0) {
                    values.add(matrix[row][col]);
                }
            }
            if (!(values.size() == 1 && values.get(0).compareTo(BigFraction.ONE) == 0)) {
                freeVariables.add(col);
            }
        }
        return freeVariables.toArray(new Integer[freeVariables.size()]);
    }
    
    public BigFraction determinant(BigFraction[][] matrix) { //method sig. takes a matrix (two dimensional array), returns determinant.
        if (matrix.length != matrix[0].length) {
            throw new RuntimeException("Must be a square matrix");
        }
        BigFraction sum = new BigFraction(0);
        int s;
        if (matrix.length == 1) {  //bottom case of recursion. size 1 matrix determinant is itself.
            return matrix[0][0];
        }
        for (int i = 0; i < matrix.length; i++) { //finds determinant using row-by-row expansion
            BigFraction[][] smaller = new BigFraction[matrix.length-1][matrix.length-1]; //creates smaller matrix- values not in same row, column
            for (int a = 1; a < matrix.length; a++) {
                for (int b = 0;b < matrix.length; b++) {
                    if (b < i) {
                        smaller[a - 1][b] = matrix[a][b];
                    }
                    else if (b > i) {
                        smaller[a - 1][b - 1] = matrix[a][b];
                    }
                }
            }
            if (i % 2 == 0) { //sign changes based on i
                s = 1;
            }
            else {
                s = -1;
            }
            sum = sum.add(new BigFraction(s).multiply(matrix[0][i].multiply(determinant(smaller)))); //recursive step: determinant of larger determined by smaller.
        }
        return sum; //returns determinant value. once stack is finished, returns final determinant.
    }
    
    public BigFraction trace(BigFraction[][] matrix) {
        if (matrix.length != matrix[0].length) {
            throw new RuntimeException("Must be a square matrix");
        }
        else {
            BigFraction sum = BigFraction.ZERO;
            for (int i = 0; i < matrix.length; i++) {
                sum = sum.add(matrix[i][i]);
            }
            return sum;
        }
    }
    
    public BigFraction[][] removeRowCol(BigFraction[][] matrix, int row, int col) {
        int n = matrix.length;
        int m = matrix[0].length;
        BigFraction[][] smaller = new BigFraction[n - 1][m - 1];
        int k = 0, l = 0;
        for (int i = 0; i < n; i++) {
            if (i == row) continue;
            for (int j = 0; j < m; j++) {
                if (j == col) continue;
                smaller[l][k] = matrix[i][j];
                k = (k + 1) % (n - 1);
                if (k == 0) l++;
            }
        }
        return smaller;
    }
    
    public BigFraction[][] matrixOfMinors(BigFraction[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length;
        BigFraction[][] matrixOfMinors = new BigFraction[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                BigFraction[][] smaller = removeRowCol(matrix, i, j);
                matrixOfMinors[i][j] = determinant(smaller);
            }
        }
        return matrixOfMinors;
    }
    
    public BigFraction[][] matrixOfCofactors(BigFraction[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length;
        BigFraction[][] matrixOfCofactors = new BigFraction[n][m];
        int sign = 1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                sign = (int)Math.pow(-1, i + j);
                matrixOfCofactors[i][j] = matrix[i][j].multiply(new BigFraction(sign));
            }
        }
        return matrixOfCofactors;
    }
    
    //returns the inverse of a matrix
    public BigFraction[][] inverse(BigFraction[][] matrix) {
        BigFraction determinant = determinant(matrix);
        if (determinant.equals(BigFraction.ZERO)) {
            System.out.println("Determinant is 0, matrix not invertible");
            return matrix;
        }
        else {
            //create matrix of minors
            BigFraction[][] matrixOfMinors = matrixOfMinors(matrix);
            //apply checkerboard to get matrix of cofactors
            BigFraction[][] matrixOfCofactors = matrixOfCofactors(matrixOfMinors);
            //transpose to get adjugate
            BigFraction[][] adjugate = transpose(matrixOfCofactors);
            //myltiply by 1/det
            BigFraction[][] inverse = multiply(adjugate, BigFraction.ONE.divide(determinant));
            return inverse;
        }
    }
    
    public BigFraction[] quadraticSolver(BigFraction a, BigFraction b, BigFraction c) {
        ArrayList<BigFraction> solutions = new ArrayList<BigFraction>();
        BigFraction two = BigFraction.TWO;
        BigFraction discriminant = b.multiply(b).subtract(new BigFraction(4).multiply(a.multiply(c)));
        if (discriminant.compareTo(BigFraction.ZERO) > 0) {
            solutions.add(((b.negate().add(new BigFraction(discriminant.pow(0.5)))).divide(two.multiply(a))));
            solutions.add(((b.negate().subtract(new BigFraction(discriminant.pow(0.5)))).divide(two.multiply(a))));
        }
        else if (discriminant.compareTo(BigFraction.ZERO) == 0) {
            solutions.add(b.negate().divide(two.multiply(a)));
        }
        else {
            System.out.println("no real solutions");
        }
        return solutions.toArray(new BigFraction[solutions.size()]);
    }
    
    public BigFraction[] eigenValues(BigFraction[][] matrix) {
        if (matrix.length != matrix[0].length) {
            throw new RuntimeException("Must be a square matrix");
        }
        else {
            BigFraction[][] identity = Matrix.identity(matrix.length);
            if (matrix.length == 2) {
                BigFraction a = BigFraction.ONE;
                BigFraction b = trace(matrix).negate();
                BigFraction c = determinant(matrix);
                //System.out.println("Characteristic polynomial of matrix is : l^2 + " + b + "l + " + c);
                return quadraticSolver(a, b, c);
            }
            else {
                System.out.println("Eigenvalue solver for matrices greater than 2x2 has not been implemented");
                return null;
            }
        }
    }
    
    public BigFraction[][] eigenVectorBasis(BigFraction[][] A, BigFraction lambda) {
        //System.out.println("A");
        //print(A);
        BigFraction[][] identity = Matrix.identity(A.length);
        BigFraction[][] IxLambda = Matrix.multiply(identity, lambda);
        //System.out.println("IxL");
        //print(IxLambda);
        BigFraction[][] IxLminusA = Matrix.subtract(IxLambda, A);
        //System.out.println("IxLminusA");
        //print(IxLminusA);
        System.out.println("basis vectors:");
        return solve(IxLminusA);
    }
    
    public BigInteger factorial(BigInteger n) {
        if (n.compareTo(BigInteger.ONE) == 0 || n.compareTo(BigInteger.ZERO) == 0) {
            return BigInteger.ONE;
        }
        else {
            return n.multiply(factorial(n.subtract(BigInteger.ONE)));
        }
    }
    
    public BigInteger combinations(BigInteger n, BigInteger k) {
        if (n.compareTo(k) < 0)
            throw new RuntimeException("n must be greater than or equal to k");
        return factorial(n).divide(factorial(k).multiply(factorial(n.subtract(k))));
    }
    
    public BigInteger permutations(BigInteger n, BigInteger k) {
        if (n.compareTo(k) < 0)
            throw new RuntimeException("n must be greater than or equal to k");
        return factorial(n).divide(factorial(n.subtract(k)));
    }
    
    public BigFraction[] expandBinomial(BigFraction a, BigFraction b, BigInteger n) {
        BigFraction[] coefficients = new BigFraction[n.add(BigInteger.ONE).intValue()];
        for (int i = 0; i < coefficients.length; i++) {
            coefficients[i] = a.pow(n.subtract(BigInteger.valueOf(i))).multiply(b.pow(i).multiply(combinations(n, BigInteger.valueOf(i))));
        }
        return coefficients;
    }
    
    public BigFraction[] square(BigFraction[] vector) {
        BigFraction[] newVector = new BigFraction[vector.length];
        for (int i = 0; i < newVector.length; i++) {
            newVector[i] = vector[i].pow(2);
        }
        return newVector;
    }
    
    public BigFraction expectedValue(BigFraction[] randomVariable) {
        BigFraction expectedValue = BigFraction.ZERO;
        for (int i = 0; i < randomVariable.length; i++) {
            BigFraction element = randomVariable[i];
            expectedValue = expectedValue.add(element.divide(randomVariable.length));
        }
        return expectedValue;
    }
    
    public BigFraction variance(BigFraction[] randomVariable) {
        return expectedValue(square(randomVariable)).subtract(expectedValue(randomVariable).pow(2));
    }
}