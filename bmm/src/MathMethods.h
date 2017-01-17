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
#include "Stack.h"
/******************************************************************************
 *                                                                            *
 *                    CONSTANTS, NUMBERS, AND VARIABLES                       *
 *                                                                            *
 ******************************************************************************/
#ifndef EPSILON
#define EPSILON 0.0000001
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef NULL
#define NULL 0
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
 
/** 
 * Function: getRow()
 * Fetch the values of the r number row and store the in the argument "store"
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

/** 
 * Function: getCol()
 * Fetch the values of the c number col and store the in the argument "store"
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

/** 
 * Function: dotProduct()
 * Compute the dot product of two vectors with length m
 * 
 * @param m the length of the vectors v1 and v2
 * 
 * @param v1 the left vector in the product
 * 
 * @param v2 the right vector in the product
 * 
 * @return the dot product of v1 and v2
 */
double dotProduct(int m, double v1[m], double v2[m]);

/** 
 * Function: euclideanNorm()
 * Compute the euclidean norm (aka magnitude or spatial length) of a given 
 * real valued m-vector
 * 
 * @param m the length of the vector v
 * 
 * @param v the vector which will be evaluated
 */
double euclideanNorm(int m, double v[m]);

/**
 *  Function: vectorProject();
 *  This function computes the projection of v1 onto v2 using the
 *  dotProduct() function. The result of the projection is stored in proj
 *
 * @param m the length of the vectors v1, v2, and proj
 *
 * @param v1 the which will be projected (projected from)
 *
 * @param v2 the vector which will recieve the projection. (projected onto)
 *
 * @param proj the stored result of the projection computation
 */
void vectorProject(int m, double v1[m], double v2[], double proj[m]);

/** 
 * Function: SwapRows()
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

/** 
 * Function: ScaleRow()
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

/** 
 * Function: elimRow()
 * This is the classic Eliminate operation. Takes in a matrix with n columns,
 * then subtracts elimBy * r2 (row 2) from r1 (row 1). We don't use this
 * function in our ref() function because of the slightly nuanced way we
 * wrote it. It works as is and I would only save a couple lines of code.
 * Not going to try to fix what works.
 * 
 * @param n the number of columns
 * 
 * @param mat the matrix to operate on
 *
 * @param r1 the row to subtract from
 * 
 * @param r2 the row to eliminate by
 * 
 * @param elimBy the scalar to multiply by r2
 * 
 */
void elimRow(int n, double mat[][n], int r1, int r2, double elimBy);

/** 
 * Function: copyMatrix()
 * Copy the values of matrix "from" into matrix "to".
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
void copyMatrix(int m, int n, double from[][n], double to[][n]);

/** 
 * Function: ref() AKA Row Echelon Form
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

/** 
 * Function: rref() AKA Reduced Row Echelon Form
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

/** 
 * Function: invert()
 * This function takes in a matrix and turns it into its multiplicative inverse.
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

/** 
 * Function: matrixMultiply()
 * This function multiplies two matrices. The left and left sides of a matrix
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
void matrixMultiply(int m, int n, int p, double left[][n], double right[][p],
        double product[][p]);

/******************************************************************************
 *                                                                            *
 *                      ANALYSIS AND CALCULUS OPERATIONS                      *
 *                                                                            *
 ******************************************************************************/

/**
 * Function: fintegral() (aka finite integral)
 * This function uses a Riemann Sum method to approximate the integral of a
 * function. 
 *
 * @param double (*funct)(double) a pointer to a function which takes a double
 *                                  and returns a double
 *
 * @param a the lower limit of integration
 *
 * @param b the upper limit of integration
 */
double fintegral(double (*funct)(double), double a, double b);

/**
 *  Function: derivative()
 *  Uses the limit definition of a derivative to approximate the 
 *  slope near a point
 *
 * @param double (*funct)(double) a pointer to a function which takes a double
 *                                  and returns a double
 *
 * @param a the point to approximate the derivative near
 */
double derivative(double (*funct)(double), double a);

//Sneak Preview: Planned functions:

//double parSum( [takes a function pointer for a sequence] );

//This one may be very hard:
//double taylor();

/******************************************************************************
 *                                                                            *
 *                              NUMERICAL METHODS                             *
 *                                                                            *
 ******************************************************************************/

/**
 * Function: binomialCoef()
 * Finds the "n choose r" value. also known as the binomial coefficient.
 * 
 * @param n the size of the set to choose from
 * @param r the number choose from a set.
 */
long int binomialCoef( long int n, long int r );

/**
 * Function: lagrangeInterpolate()
 * This function will make a polynomial of degree n-1 given n points. It uses
 * an algorithm known as lagrange interpolation. The function executes in the 
 * following manner:
 *
 *        (let n-1 be the order/degree for simplicity) 
 *              
 *        1) -> Take n pairs of  x and y data points
 *        2) -> For the set of points, load them into a matrix like so:
 *              
 *              [  (1)    (x[0])    (x[0])^2   ...   (x[0])^(d)  ]
 *              [  (1)    (x[1])    (x[1])^2   ...   (x[1])^(d)  ]
 *         X =  [                               .                ]
 *              [                               .                ]
 *              [                               .                ]
 *              [  (1)   (x[n-1])  (x[n-1])^2  ...  (x[n-1])^(d) ]
 *
 *
 *        3) -> Use the array of y data as a column vector.
 *
 *                              [  y[0]  ]
 *                              [  y[1]  ]
 *                         y =  [    .   ]
 *                              [    .   ]
 *                              [    .   ]
 *                              [ y[n-1] ]
 *
 *
 *        4) -> Use the array c as the vector of unknown coefficients which
 *              we want to solve for. We set up the matrix equation:
 *
 *                                 Xc = y
 *
 *
 *        5) -> Since X is invertible by design, it has a solution for c :
 *
 *                              c = (X^(-1))y
 *
 *
 * @param n the number of data points. (one higher than the degree)
 *
 * @param x[n] the set of x data points
 *
 * @param y[n] the set of y data points
 *
 * @param c[n] the array which will store all the coefficients of the
 *               the polynomial. c[0] is the constant term, c[d] is the highest
 *               power's coefficient.
 */
void lagrangeInterpolate(long int n, double x[n], double y[n], double c[n]);

/**
 * Function: meanInterpolate()
 * This function will make a best-fit polynomial of a specified degree. Given 2
 * vectors representing the x and y coordinates of a graph, the function will
 * analyze them in the following manner:
 *             
 *        1) ->  Call the interpolate function using the given x and y vectors
 *
 *        2) -> Do this for recursively for every n size subset of our data 
 *              points, taking the mean of the results as you compute. The 
 *              result is a vector of coefficients for a best fit polynomial of 
 *              degree d.
 *  
 *  (NOTE)  2)-> Recursive computation finds a num_coefs size subset of the x,y
 *               data pairs then calls lagrangeInterpolate() on them to compute.
 *               only guaranteed to work if you pass -1 to path.
 *
 * @param path an iterator for recursive computation. 
 *                ALWAYS PASS AS -1 TO INITIALIZE!
 * 
 * @param *combo this is for combinatorially iterating through possible 
 *                  indeces when interpolating. 
 *                  ALWAYS PASS AS NULL POINTER TO INITIALIZE!
 *
 * @param num_points the number of (x, y) data points.
 * 
 * @param x[num_points] the set of x data points
 *
 * @param y[num_points] the set of y data points
 *
 * @param num_coefs one more than the degree of a polynomial you want
 *
 * @param c[num_coefs] the array which will store all the coefficients of the
 *               the polynomial. c[0] is the constant term, c[d] is the highest
 *               power's coefficient.
 */
void meanInterpolate(long int path, Stack *combo, long int num_points, 
                    double x[num_points], double y[num_points], 
                    long int num_coefs, double c[num_coefs]);
#endif /* MATHMETHODS_H */
