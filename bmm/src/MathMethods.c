#include "../hed/MathMethods.h"
#include <stdlib.h>
#include <math.h>

/******************************************************************************
 *                                                                            *
 *                    CONSTANTS, NUMBERS, AND VARIABLES                       *
 *                                                                            *
 ******************************************************************************/

//unless otherwise defined, delta is a small number near zero
static double delta = EPSILON;

/*This function is used to set a static system variable called delta.
 *delta is used to take integrals, derivatives, and other values
 */
void setDelta(double d)
{
    delta = d;
}

/******************************************************************************
 *                                                                            *
 *                         LINEAR ALGEBRA OPERATIONS                          *
 *                                                                            *
 ******************************************************************************/

/* Fetch the values of the r number row and store the in the argument "store"
 * 
 * @param m the number of rows of the matrix mat
 * 
 * @param n the number of cols of the matrix mat
 * 
 * @param mat the matrix you wish to read a row from
 * 
 * @param r the row you want to read
 * 
 * @the store the vector you wish to store into
 */
void getRow(int m, int n, double mat[][n], int r, double store[n])
{
    int j;
    for (j = 0; j < n; j++) {
        store[j] = mat[r][j];
    }
    return;
}

/* Fetch the values of the c number col and store the in the argument "store"
 * 
 * @param m the number of rows of the matrix mat
 * 
 * @param n the number of cols of the matrix mat
 * 
 * @param mat the matrix you wish to read a column from
 * 
 * @param c the column you want to read
 * 
 * @the store the vector you wish to store into
 */
void getCol(int m, int n, double mat[][n], int c, double store[n])
{
    int i;
    for (i = 0; i < m; i++) {
        store[i] = mat[i][c];
    }
    return;
}

/* Compute the dot product of two vectors with length m
 * 
 * @param m the length of the vectors v1 and v2
 * 
 * @param v1 the left vector in the product
 * 
 * @param v2 the right vector in the product
 * 
 * @return the dot product of v1 and v2
 */
double dotProduct(int m, double v1[m], double v2[m])
{
    double dot = 0;
    int i;
    for (i = 0; i < m; i++) {
        dot += v1[i] * v2[i];
    }
    return dot;
}

/* Function: SwapRows()
 * This function takes in an n-column matrix. We assume that row1 and row2 are 
 * valid rows inside of your n-column matrix. We swap the values within row1 and
 * row2. This means the matrix gets modified.
 *
 * @param n number of cols of your matrix
 * 
 * @param mat your matrix pointer
 * 
 * @param row1 your first row
 * 
 * @param row2 your other row
 */
void swapRows(int n, double mat[][n], int r1, int r2)
{
    int j;
    double temp;
    for (j = 0; j < n; j++) {
        temp = mat[r1][j];
        mat[r1][j] = mat[r2][j];
        mat[r2][j] = temp;
    }
    return;
}

/* Function: ScaleRow
 * Take in an n-column matrix and scale the the r row of your matrix by a double
 * labeled scaleBy.
 * 
 * @param n the number of columns
 * 
 * @param mat your matrix pointer
 * 
 * @param r the row you want to modify
 * 
 * @param scaleBy the scalar you wish to multiply the row by
 * 
 */
void scaleRow(int n, double mat[][n], int r, double scaleBy)
{
    int j;
    for (j = 0; j < n; j++) {
        mat[r][j] *= scaleBy;
    }
    return;
}

/* This is the classic Eliminate operation. Takes in a matrix with n columns,
 * then subtracts elimBy * r2 (row 2) from r1 (row 1). We don't use this
 * function in our ref() function because of the slightly nuanced way we
 * wrote it. It works as is and I would only save a couple lines of code.
 * Not going to try to fix what works.
 * 
 * @param n the number of columns
 * 
 * @param r1 the row to subtract from
 * 
 * @param r2 the row to eliminate by
 * 
 * @param elimBy the scalar to multiply by r2
 * 
 */
void elimRow(int n, double mat[][n], int r1, int r2, double elimBy)
{
    int j;
    for (j = 0; j < n; j++) {
        mat[r1][j] -= elimBy * mat[r2][j];
    }
    return;
}

/* Copy the values of matrix "from" into matrix "to".
 * 
 * @param m the number of rows
 * 
 * @param n the number of cols
 * 
 * @param from the matrix with desired values
 * 
 * @param to the matrix to copy into
 *
 */
void copy(int m, int n, double from[][n], double to[][n])
{
    int i;
    int j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            to[i][j] = from[i][j];
        }
    }
}

/* Function: ref() AKA Row Echelon Form
 * This is an implementation of the well-known Gaussian Elimination
 * algorithm taught in an introductory Linear Algebra course. It takes in a 
 * matrix of "doubles," which are high-precision finite decimal expansions of 
 * Real numbers. The next two arguments to the function are the number of rows, 
 * m, and the number of columns, n, of the matrix of Real numbers. The prefix 
 * "unsigned int" means that I am restricting the inputs to the integer 
 * numbers, as defined in C. The function returns the determinant if it was 
 * square and zero otherwise. We assume you know what kind of matrix you are 
 * passing it.
 * 
 * @param m the number of rows of mat
 *
 * @param n the number of columns of mat
 * 
 * @param mat a pointer to the upper leftmost address (in physical computer 
 *          memory of a matrix. mat will be modified by this algorithm, so only 
 *          give this function a matrix you want to change.
 * 
 * @return if the matrix was square, return the determinant.
 *          if the matrix was not square, we return a zero, which is just a 
 *          placeholder value so that we dont return garbage.
 */
double ref(int m, int n, double mat[][n])
{
    int x = 0;
    int y = 0;
    //these two numbers represent the x row and y column of the matrix
    //these are our iterators
    int i = 0;
    //this is the top row of the jth column. otherwise known as the ith row.
    int j = 0;
    //this is going to be the first nonzero column, starting from the left
    //when this flag goes up, our matrix did 
    int noPivotColFlag = 0;
    //this next int flag is for when our whole row is a zero row
    double det = 1;
    //check if our matrix is square.
    if (m != n) {
        det = 0;
    }
    while (i < m && j < n) {
        /* due to floating point rounding, we may end up with the difference
         * between a number and another being nonzero, even when by our
         * calculations by hand, they are equal. The constant
         * EPSILON is a buffer to make sure we still treat zeros
         * properly. A note: fabs(x) is a function such that it takes a real
         * number and returns the real number's magnitude. "fabs" is short for
         * "floating point absolute value"
         */
        if (fabs(mat[i][j]) < EPSILON) {//if matrix's ij value is near zero
            for (x = i; x < m; x++) {//iterate through rows (vertically)
                if (fabs(mat[x][j]) > EPSILON) {
                    //if we find a nonzero entry below, swap them
                    swapRows(n, mat, i, x);
                    det *= -1; //det -> negative det after a swap
                    break;
                } else if (x == m - 1) {
                    //if we reach the bottom and they were all zero
                    //then the column we were looking at is all zeros
                    //below row i. We turn on a variable to signify this
                    noPivotColFlag = 1;
                }
            }
        }
        //this checks if we just iterated through a column 
        //without a pivot position
        if (noPivotColFlag == 1) {
            noPivotColFlag = 0; //if we did, turn the flag off
            j++; //go to the next column
            det = 0; //we now know the determinant is zero
            continue;
        }
        //we scale our current row next. adjust the determinant accordingly
        det *= mat[i][j];
        scaleRow(n, mat, i, 1 / mat[i][j]);
        //we now do an elimination operation to each of the rows below.
        //the upper limits of these loops
        for (x = i + 1; x < m; x++) {
            for (y = n - 1; y >= j; y--) {
                //due to order of computation, this loop is computer backwards.
                //this is a rather advanced C topic. Think about it for a bit
                //and perhaps it makes a bit more sense.
                mat[x][y] -= mat[i][y] * mat[x][j];
            }
        }
        j++;
        i++; //all good, iterate forward!
    }
    //all done, return the determinant
    return det;
}

/* Function: rref() AKA Reduced Row Echelon Form
 * This function calls ref() once to make sure that it toying with a matrix
 * in row echelon form. After that it starts iterating from row (m - 1) all 
 * the way to the top, eliminating the values of the leading ones as it goes.
 * 
 * Keep in mind that if you write:
 * 
 *      double mat[m][n] = {...some set of vectors...};
 *      ref(m, n, mat);
 *      double det = rref(mat);
 * 
 * Then double now stores 1. In fact, if you call ref() or rref() on any matrix
 * more than once, you make that matrix have a determinant of 1.
 * 
 * @param m the number of rows
 * 
 * @param n the number of columns
 * 
 * @param mat the matrix which you wish to eliminate on. 
 * 
 * @return If the matrix was square, return the determinant. if it was not,
 *          return zero. We assume you understand that this zero has no
 *          correlation to a determinant for the m not equal to n case.
 */
double rref(int m, int n, double mat[][n])
{
    double det = ref(m, n, mat);
    //these are our loop iterators, x and y:
    int x;
    int y;
    //these are the values of the row/column we are working on
    int i = m - 1;
    int j;
    //first find the first nonzero row from the bottom up:
    //this is the thing that tells us we found a pivot column
    int foundPivotFlag = 0;
    while (i >= 0) {
        //first we find the first pivot position, from the bottom up
        for (x = i; x >= 0; x--) {
            for (y = 0; y < n; y++) {
                if (mat[x][y] > EPSILON) {
                    i = x;
                    j = y;
                    foundPivotFlag = 1;
                    break;
                }
            }
            if (foundPivotFlag == 1) {
                foundPivotFlag = 0;
                break;
            }
        }
        //[i][j] is the location of the bottom most pivot position now
        //Now we have to do an elimination operation, going up the matrix
        for (x = i - 1; x >= 0; x--) {
            for (y = n - 1; y >= j; y--) {
                //due to order of computation, this loop is computer backwards.
                //this is a rather advanced C topic. Think about it for a bit
                //and perhaps it makes a bit more sense.
                mat[x][y] -= mat[i][y] * mat[x][j];
            }
        }
        i--;
    }
    //All Done, return the determinant
    return det;
}

/* This function takes in a matrix and turns it into its multiplicative inverse.
 * It copies the matrix into an augmented matrix which is twice the width of the 
 * original. The function then does an rref() on the augmented matrix. After 
 * that it takes the right half of the reduced augmented matrix and copies it
 * the matrix supplied by the argument.
 * 
 * @param m the number of rows
 * 
 * @param n the number of columns
 * 
 * @param return det(mat) if the matrix was invertible, return 0 otherwise
 * 
 */
double invert(int m, int n, double mat[][n])
{
    double det = 0;
    int i;
    int j;

    if (m != n) {
        return 0;
    }

    double checkMat[m][n];
    //use this matrix to check what the determinant is before we try to find
    //and inverse.

    copy(m, n, mat, checkMat);

    det = ref(m, n, checkMat);

    if (fabs(det) < EPSILON) {
        return 0;
    }

    double augMat[m][2 * n];

    //instead of using the copy function, we must use a custom loop
    //this is because copy function requires 2 matrices of the same dimensions

    for (i = 0; i < m; i++) {
        for (j = 0; j < 2 * n; j++) {
            if (j >= n) {
                augMat[i][j] = 0;
                if ((j - n) == i) { //the we are in the augmented diagonal
                    augMat[i][j] = 1;
                }
            } else {

                augMat[i][j] = mat[i][j];
            }
        }
    }

    rref(m, n * 2, augMat);
    for (i = 0; i < m; i++) {
        for (j = n; j < 2 * n; j++) {
            mat[i][j - n] = augMat[i][j];
        }
    }

    return det;

}

/* This function multiplies two matrices. The left and left sides of a matrix
 * product must share a corresponding column/row match up. That is to say that
 * the columns of the left side of the product are the same as the rows of the
 * right side of the product. This is represented by the argument n. The result
 * of the matrix multiplication is stored in the matrix labeled "product."
 * 
 * @param m the rows of the left matrix in the product
 * 
 * @param n the columns of the left and the rows of the right in a matrix
 * 
 * @param p the columns of the right in the product
 * 
 * @param left the left side of the product
 * 
 * @param right the right side of the product
 * 
 * @param product the result of our multiplication
 * 
 */
void multiply(int m, int n, int p, double left[][n], double right[][p],
        double product[][p])
{
    int i;
    int j;
    double *v1 = calloc(n, sizeof(double));
    double *v2 = calloc(n, sizeof(double));;
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            getRow(m, n, left, i, v1);
            getCol(n, p, right, j, v2);
            product[i][j] = dotProduct(n, v1, v2);
        }
    }
    free(v1);
    free(v2);
}

/******************************************************************************
 *                                                                            *
 *                      ANALYSIS AND CALCULUS OPERATIONS                      *
 *                                                                            *
 ******************************************************************************/

/*
 *finite integral
 */
double fintegral(double (*funct)(double), double a, double b)
{
    //value of a reimann sum
    double riemann = 0;
    //index to calulate an integral with
    double index;
    //a flag and swap value for when you integrate from higher to lower
    double backwards = 1;

    //use negative of an integral property
    if (b < a) {
        backwards = b;
        b = a;
        a = backwards;
        backwards = -1;
    }

    //start indexing
    index = a;

    //calulate the riemann sum
    while (index < b) {
        riemann += funct(index) * delta;
        index += delta;
    }

    //multiply by the factor mentioned to get the finite integral approximation
    return riemann * backwards;
}

/*
 * approximate derivative
 */
double derivative(double (*funct)(double), double a)
{
    return (funct(a + delta) - funct(a)) / delta;
}

/******************************************************************************
 *                                                                            *
 *                              NUMERICAL METHODS                             *
 *                                                                            *
 ******************************************************************************/

/*
 * This function will make a fit polynomial of a specified degree. Given 2
 * vectors representing the x and y coordinates of a graph, the function will
 * analyze them in the following manner:
 *             
 *        (let d be the degree for simplicity) 
 *              
 *        1) -> Take (degree + 1) pairs of  x and y data points
 *        2) -> Pick d
 *              Load the x[j] data into a matrix like so: 
 *              
 *              [  (1)    (x[1])    (x[1])^2   ...   (x[1])^(d)  ]
 *              [  (1)    (x[2])    (x[2])^2   ...   (x[2])^(d)  ]
 *         X =  [                               .                ]
 *              [                               .                ]
 *              [                               .                ]
 *              [  (1)   (x[d+1])  (x[d+1])^2  ...  (x[d+1])^(d) ]
 *
 *
 *        3) -> Use the array of y[j] data as a column vector.
 *
 *                              [  y[1]  ]
 *                              [  y[2]  ]
 *                         y =  [    .   ]
 *                              [    .   ]
 *                              [    .   ]
 *                              [ y[d+1] ]
 *
 *        4) -> Use the array c as the vector of unknown coefficients which
 *              we want to solve for. We set up the matrix equation:
 *
 *                                 Xc = y
 *
 *        5) -> Since X is invertible by design, it has a solution for c:
 *
 *                              c = (X^(-1))y
 *
 *        6) -> We do this for the first d+1 points, the second d+1 points, 
                until eventually we have n - d polynomials.
 *              After that, take the mean of all of the results. The result
 *              is a vector of coefficients for a fit polynomial of 
 *              degree d.
 *
 * @param n the number of (x, y) data points.
 * 
 * @param x[n] the set of x data points
 *
 * @param y[n] the set of y data points
 *
 * @param d the degree of a polynomial you want
 *
 * @param c[d+1] the array which will store all the coefficients of the
 *               the polynomial. c[0] is the constant term, c[d] is the highest
 *               power's coefficient. This will be filled with all zeroes if
 *               d + 1 > n.
 * @return void
 */
void interpolate(int n, double x[n], double y[n],
                 int d, double c[d + 1])
{   
    //loop iterators

    int i = 0;
    int j = 0;
    int k = 0;
    
    //this is the number of coefficients. it is also the length of c
    //such a variable is much easier to read than "d+1" over and over.
    const int coef = d+1; 
    //matrix of x values raised to powers
    double xMat[coef][coef];
    //intialize xMat
    for(i = 0; i < coef; i++)
    {
        for(j = 0; j < coef; j++)
        {
            //printf("\nyou got this far %d\n", j);
            xMat[i][j] = 0;
        } 
    }

    //column vector of y values
    double yMat[coef][1];
    //initialize yMat
    
    for(j = 0; j < coef; j++)
    {
        yMat[coef][0] = 0;
    }
    
    //column vector which is the product of inverse(xMat) and yMat   
    double pMat[coef][1];
    //intialize xMat
    for(i = 0; i < coef; i++)
    {
        pMat[i][0] = 0;
    }
   
    // check if we have enough data to build the polynomial
    if(d+1 > n){
        //if we don't that's a no-no. set c to zero vector and exit.
        for(i = 0; i < coef; i++)
        {
            c[i] = 0;
        }
        return;
    }
    
    //start building the polynomials. there are n - d of them.
    for(i = 0; i < n - d; i++)
    {   
        //build a single polynomial
        for(j = 0; j < d+1 + 0; j++)
        {
            yMat[j][0] = y[j + i];
            //printf(" %+2.4f ", yMat[j][0]);
            for(k = 0; k < coef; k++)
            {
                xMat[j][k] = x[j + i];
                xMat[j][k] = powl(xMat[j][k], k);
                
                //printf(" %+2.4f ", xMat[j][k]);
            }
            //printf("\n");
              
        }
        //take the inverse of the xMat matrix so we can find the 
        //coefficients of our polynomial

        
        invert(coef, coef, xMat);
        //multiply this inverse matrix by yMat to get  one of the 
        //polynomial coefficient vectors
        
        multiply(coef, coef, 1, xMat, yMat, pMat);

        //c has i previous vectors stored in it through an average.
        //in order for it to have the "weight" of i iterations, we
        //multiply it by i
        
        for(k = 0; k < coef; k++)
        {
                c[k] = (i*c[k] + pMat[k][0]) / (i + 1); 
                //denominator is i + 1 because of the weighted average
        }
    }
    //done
    return;
}
