/* 
 * File:   MathMethods.h
 * Author: Abram Rodgers
 * 
 * This is a computational methods library. Somewhat meant to exhibit the kinds
 * of math I am familiar with, but also for useful application in the future.
 * 
 *
 * Created on March 10, 2016, 12:07 AM
 */

#ifndef MATHMETHODS_H
#define MATHMETHODS_H

#include <stdint.h>

/******************************************************************************
 *                                                                            *
 *                    CONSTANTS, NUMBERS, AND VARIABLES                       *
 *                                                                            *
 ******************************************************************************/
#ifndef EPSILON
#define EPSILON 0.000001
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*This function is used to set a static system variable called delta.
 *delta is used to take integrals, derivatives, and other values
 */
void setDelta(double d);

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
void getRow(int m, int n, double mat[][n], int r, double store[n]);

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
void getCol(int m, int n, double mat[][n], int c, double store[n]);

//write the documentation!
double dotProduct(int m, double v1[m], double v2[m]);

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
void swapRows(int n, double mat[][n], int r1, int r2);

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
void scaleRow(int n, double mat[][n], int r, double scaleBy);


/* This is the classic Eliminate operation. Takes in a matrix with n columns,
 * then subtracts elimBy * r2 (row 2) from r1 (row 1).
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
void elimRow(int n, double mat[][n], int r1, int r2, double elimBy);

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
void copy(int m, int n, double from[][n], double to[][n]);

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
double ref(int m, int n, double mat[][n]);

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
double rref(int m, int n, double mat[][n]);

/* This function takes in a matrix and turns it into its multiplicative inverse
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
double invert(int m, int n, double mat[][n]);

/* This function multiplies two matrices. The left and left sides of a matrix
 * product must share a corresponding column/row match up. That is to say that
 * the columns of the left side of the product are the same as the rows of the
 * right side of the product. This is represented by the argument n. The result
 * of the matrix multiplication is stored in the matrix labeled "product."
 * 
 * @param m the rows of the left matrix in a product
 * 
 * @param n the columns of the left and the rows of the right in a matrix
 * 
 * @param p the columns of the right of the product
 * 
 * @param left the left side of the product
 * 
 * @param right the right side of the product
 * 
 * @param product the result of our multiplication
 */
void multiply(int m, int n, int p, double left[][n], double right[][p],
        double product[][p]);

/******************************************************************************
 *                                                                            *
 *                      ANALYSIS AND CALCULUS OPERATIONS                      *
 *                                                                            *
 ******************************************************************************/

/* The following functions need documentation written for them
 */

//fintegral? finite integral. it is a pun. (Das war ein witz.)
//documentation to be written
double fintegral(double (*funct)(double), double a, double b);

//approximate slope near a point
//documentation to be written
double derivative(double (*funct)(), double a);

//double parSum( [takes a function pointer for a sequence] );

//This one may be very hard:
//double taylor();

/******************************************************************************
 *                                                                            *
 *                              NUMERICAL METHODS                             *
 *                                                                            *
 ******************************************************************************/

/*
 * This function will make a best-fit polynomial of a specified degree. Given 2
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
 *        6) -> Do this for every d+1 size subset of our data points.
 *              After that, take the mean of all of the results. The result
 *              is a vector of coefficients for a best fit polynomial of 
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
 *               power's coefficient.
 */
void interpolate(int n, double x[n], double y[n],
                 int d, double c[d + 1]);
#endif /* MATHMETHODS_H */
