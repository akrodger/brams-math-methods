/* 
 * File:   MathMethods.h
 * Author: Abram Rodgers
 * Email: aks.rodgers@gmail.com
 * 
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
 *                       MEMORY MANAGEMENT OF MATRICES                        *
 *                                                                            *
 ******************************************************************************/

/**
 * Macro: initMat()
 * just a shortcut for doing calloc() on a matrix. Calloc initializes values to
 * zero. 
 *
 * A Word on how memory is managed with matrices and vectors:
 * A m and n matrix is a n m by n 2D array of doubles with C99 representation.
 * An m by 1 matrix is a column vector.
 * Let vec be a column vector. the vec[0] is an an array of length n, almost.
 * It is the memory location of the the column vector, and therefore may be
 * used as such. An exploitation of this fact is in the matrixMultiply function.
 *
 * Actual definition is in header.
 * 
 * @param m the number of rows of the matrix mat
 * 
 * @param n the number of cols of the matrix mat
 * 
 * @param mat the matrix you wish to initialize
 */
#ifndef initMat
#define initMat(m, n, mat) mat = calloc((m)*(n), sizeof(double))
#endif

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
 * @param store the array you wish to store into
 */
void getRow(int m, int n, double (*mat)[n], int r, double store[n]);

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
 * @the store the array you wish to store into
 */
void getCol(int m, int n, double (*mat)[n], int c, double store[n]);

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
double dotProduct(int m, double (*v1)[1], double (*v2)[1]);

/** 
 * Function: euclideanNorm()
 * Compute the euclidean norm (aka magnitude or spatial length) of a given 
 * real valued m-vector
 * 
 * @param m the length of the vector v
 * 
 * @param v the vector which will be evaluated
 */
double euclideanNorm(int m, double (*v)[1]);

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
void vectorProject(int m, double (*v1)[1], double (*v2)[1], double (*proj)[1]);

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
void swapRows(int n, double (*mat)[n], int r1, int r2);

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
void scaleRow(int n, double (*mat)[n], int r, double scaleBy);

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
void elimRow(int n, double (*mat)[n], int r1, int r2, double elimBy);

/** 
 * Function: copyMatrix()
 * Copy the values of matrix "from" into matrix "to".
 * You can use this as a reguar array copy by setting n=1 and passing 2 arrays
 * C might get mad at you though and send all those "warnings"
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
void copyMatrix(int m, int n, double (*from)[n], double (*to)[n]);

/**
 * Function: transpose()
 * The transpose of an m by n matrix A is the n by m matrix A' such that the
 * elements of the rows of A are the elements of the columns of A'.
 * 
 * Other single-matrix manipulation functions in this library directly
 * manipulate the given matrix. (See ref and invert functions as examples.)
 * However an m by n and an n by m 2D array might not be handled the same way
 * in memory, depending on what processor one uses. (E.g. I compiled this
 * library on a microcontroller and weird memory stuff happened with index out
 * of bounds.) For this reason, we require that the matrix we wish to transpose
 * and the storage location of the transposed matrix are different.
 *
 * However, it is worth noting that memory allocation for m by n and n by m are
 * identical. In theoy, it could work so that transpose(m, n, A, A); is valid.
 * I don't think this is worth implementing though, as it would need a swapping
 * mechanism. If you disagree email me, it'd interesting to know what you think.
 *
 * @param m the number of rows in mat.
 * 
 * @param n the number of columns mat.
 *
 * @mat the matrix to transpose.
 *
 * @matTranspose the matrix to store the transpose in.
 */
void transpose(int m, int n, double (*mat)[n], double (*matTranspose)[m]);

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
double ref(int m, int n, double (*mat)[n]);

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
double rref(int m, int n, double (*mat)[n]);

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
double invert(int m, int n, double (*mat)[n]);

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
        double (*product)[p]);
/**
 * Function: arrayAdd()
 * Add the corresponding elements of two 1D arrays of doubles. 
 * Return value is the last argument.
 *
 * @param m the length of the arrays.
 *
 * @param left the left addend
 *
 * @param right the right addend
 *
 * @param sum the sum of the left and the right
 */
void arrayAdd(int m, double *left, double *right, double *sum);
/**
 * Function: arrayScale()
 * Multiplies all elements of an array of doubles by a number.
 *
 * @param m the length of the arrays.
 *
 * @param arr the array to be scaled
 *
 * @param scaleBy the double with which you scale the vector
 *
 * @param scaled the scaled array
 */
void arrayScale(int m, double *arr, double scaleBy, double *scaled);
/**
 * Function: arrayAdd()
 * Add the corresponding elements of two matrices. 
 * Return value is the last argument.
 *
 * @param m the number of rows.
 *
 * @param n the number of columns
 *
 * @param left the left addend
 *
 * @param right the right addend
 *
 * @param sum the sum of the left and the right
 */
void matrixAdd(int m, int n, double (*left)[n], double (*right)[n], 
				double (*sum)[n]);
/******************************************************************************
 *                                                                            *
 *                 NUMERICAL ANALYSIS AND CALCULUS OPERATIONS                 *
 *                                                                            *
 ******************************************************************************/

/**
 * Function: rectIntegral() 
 * This function uses a Riemann Sum method with fixed width size over rectangles
 * to approximate the integral of a function mapping from doubles to doubles. 
 *
 * @param double (*funct)(double) a pointer to a function which takes a double
 *								  and returns a double
 *
 * @param a the lower limit of integration
 *
 * @param b the upper limit of integration
 */
double rectIntegral(double (*funct)(double), double a, double b);

/**
 * Function: simpIntegral() 
 * This function uses integrated interpolating parabolas
 * to approximate the integral of a function mapping from doubles to doubles.
 * the static double DELTA is used as the assumed separation between integrating
 * parabolas.
 *
 * Technically, this is the Cavalieri-Simpson rule, meaning we use averaged
 * midpoints. 
 *
 * Formula named for Thomas Simpson (1710-1761) but Kepler was using similar
 * methods 100 years earlier. (Citation: Wikipedia page on Simpson Rule)
 *
 * @param double (*funct)(double) a pointer to a function which takes a double
 *								  and returns a double
 *
 * @param a the lower limit of integration
 *
 * @param b the upper limit of integration
 */
double simpIntegral(double (*funct)(double), double a, double b);

/**
 *  Function: simpleDerivative()
 *  Uses the limit definition of a derivative to approximate the 
 *  slope near a point. Literally only 1 line as a 1st order approx.
 *
 * @param double (*funct)(double) a pointer to a function which takes a double
 *								  and returns a double
 *
 * @param a the point to approximate the derivative near
 */
double simpleDerivative(double (*funct)(double), double a);

//Sneak Preview: Planned functions:

//double parSum( [takes a function pointer for a sequence] );

//This one may be a little hard:
//double taylor();

/**
 * This function uses Adams Bashforth 3-step method to solving systems of 
 * systems of ordinary differential equations. A derivation for AB3 is on
 * wikipedia. Points are saved to the argument y_initial.
 *
 * Words of Caution: 2 things. First. numPnts = total number of solves + 1. This
 *					 is because we don't need to solve for initial point. Also,
 *					 note that in the case of (iostep == numPnts) being true, we
 *					 actually simply return the last point in the iteration. So 
 *					 if Itell the solver to compute 5000 points and saved to
 *					 y_solve once every 5000 poinths, then y_solve will be an n 
 *					 by 1 column vector containing the last point computed.
 *
 * @param n The dimension of the differential equation.
 * 
 * @param void (*vecField)(double[n], double[n])) Function pointer to the
 * 				vector/flow field. The first argument is the input to the
 *				field, the second argument is the output.
 *
 * @param (*y_initial)[1] an n-dimensional column vector which gives the 
 * 							initial condition of the ODE
 *
 * @param deltaT The time step for the solver to use
 * 
 * @param iostep The ode saves a point to the solution once every iostep
 *					number of iterations. Make sure to make this larger if
 *					you want to use a really small deltaT.
 *
 * @param numPnts This is the number of points the solver finds. Must be at
 * 					at least 
 *
 * @param (*y_solved)[(numPnts - (numPnts % iostep))/iostep] 
 *						The matrix containing the ODE numerical
 *						solution. Every column is a point on the trajectory.
 *
 * @param *time an array containing the time at each point.
 */
void adamsBash3(int n, void (*vecField)(double[n], double[n]),
				double (*y_initial)[1], double deltaT, 
				int iostep, long int numPnts,
				double (*y_solve)[(numPnts - (numPnts % iostep))/iostep], 
				double *time);

/**
 * Function: binomialCoef()
 * Finds the "n choose r" value. also known as the binomial coefficient.
 * 
 * @param n the size of the set to choose from
 * @param r the number choose from a set.
 */
long int binomialCoef( long int n, long int r );

/**
 * Function: vanderInterpolate()
 * This function will make a polynomial of degree n-1 given n points. It uses
 * an algorithm known as Vandermonde Interpolation. The function executes in the 
 * following manner:
 *
 *		(let n-1 be the order/degree for simplicity) 
 *			  
 *		1) -> Take n pairs of  x and y data points
 *		2) -> For the set of points, load them into a matrix like so:
 *			  
 *              [  (1)    (x[0])    (x[0])^2   ...   (x[0])^(d)  ]
 *              [  (1)    (x[1])    (x[1])^2   ...   (x[1])^(d)  ]
 *         X =  [                               .                ]
 *              [                               .                ]
 *              [                               .                ]
 *              [  (1)   (x[n-1])  (x[n-1])^2  ...  (x[n-1])^(d) ]
 *
 *
 *		3) -> Use the array of y data as a column vector.
 *
 *                              [  y[0]  ]
 *                              [  y[1]  ]
 *                         y =  [    .   ]
 *                              [    .   ]
 *                              [    .   ]
 *                              [ y[n-1] ]
 *
 *
 *		4) -> Use the array c as the vector of unknown coefficients which
 *			  we want to solve for. We set up the matrix equation:
 *
 *								 Xc = y
 *
 *
 *		5) -> Since X is invertible by design, it has a solution for c :
 *
 *							  c = (X^(-1))y
 *
 *
 * @param n the number of data points. (one higher than the degree)
 *
 * @param x[n] the set of x data points
 *
 * @param y[n] the set of y data points
 *
 * @param c[n] the array which will store all the coefficients of the
 *			   the polynomial. c[0] is the constant term, c[d] is the highest
 *			   power's coefficient.
 */
void vanderInterpolate(long int n, double x[n], double y[n], double c[n]);

/**
 * Function: lagrangeInterpolate()
 * This function calculates the same polynomial which the Vandermonde method
 * does, but instead it uses the less numerically touchy method of computing.
 * This method uses a set of lagrange polynomials as a basis for the
 * n-dimensional double-valued polynomial vector space. The computation
 * does not involve inverting a matrix, where as the Vandermonde method does.
 *
 * For a set of n data points, we construct a polynomial of degree n-1
 * where the polynomial is a linear combination of a lagrange basis formed 
 * with the given set of data points.
 * 
 * The kth memeber of the basis l[k] is 
 *
 *		l[k](s) = PODUCT(j=0 to j=n-1, skip j=k) (s - x[j])/(x[k] - x[j])
 *
 * The final polynomial F(s) is a linear combination of this basis:
 * 
 * 		F(s) = SUM(k=0 to k=n-1) y[k] * l[k](s)
 *
 * @param n the number of data points. (one higher than the degree)
 *
 * @param x[n] the set of x data points
 *
 * @param y[n] the set of y data points
 *
 * @param c[n] the array which will store all the coefficients of the
 *			   the polynomial. c[0] is the constant term, c[d] is the highest
 *			   power's coefficient.
 */
//Not Yet Implemented
//void lagrangeInterpolate(long int n, double x[n], double y[n], double c[n]);

/**
 * Function: meanPolynomial()
 * This function will make a "best-fit" polynomial of a specified degree. Given
 * 2 vectors representing the x and y coordinates of a graph, the function will
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
 *               data pairs then calls vanderInterpolate() on them to compute.
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
void meanPolynomial(long int path, Stack *combo, long int num_points, 
                    double x[num_points], double y[num_points], 
                    long int num_coefs, double c[num_coefs]);

/**
 * Function: discreteLeastSquares()
 * The "meanPolynomial" function above computes the average polynomial of 
 * a certain degree with repect to a data set. Here's a more general algorithm.
 *
 * Where as above, we had a polynomial of degree num_coefs-1. For brevity, say
 * num_coefs-1 = m. Therefore, we can consider the polynomial to be a set of
 * m functions. The goal of Least Squares approximation is to represent a data
 * set with a sum of many functions. Moreover, we don't want to iterate through
 * a huge number of subsets to do it. (Like in meanPolynomial.) We use linear
 * algebra to solve this problem. 
 *
 * Core idea:
 * double (*phi[m])(double) is an array of m mappings from doubles to doubles.
 * Assume they are linearly inderpendent of one another. Now sample each of the 
 * n many x-axis data points at every function. Store these in an m by n matrix.
 * Call this matrix B. Let B' be the transpose of B. Assume y is a vector with
 * all the y data points. The coordinate vector of the data set with respect to
 * the phi functions will be:
 *								a = (((B'B)^-1)B')y
 * The proof of this fact requires some rather advanced applications of linear
 * algebra an some multivariable calculus. I did my undergraduate thesis on this
 * subject, you can email me for my writing on this. It is also on Wikipedia,
 * although the Wikipedia article is a little dense.
 */
void discreteLeastSquares(int n, double x[n], double y[n], int m,
							double (*phi[m])(double), double result[m]);
#endif /* MATHMETHODS_H */
