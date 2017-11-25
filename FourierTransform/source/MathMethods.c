/* 
 * File:   MathMethods.c
 * Author: Abram Rodgers
 * Email@ aks.rodgers@gmail.com 
 *
 * This is a computational methods library. Somewhat meant to exhibit the kinds
 * of math I am familiar with, but also for useful application in the future.
 * 
 *
 * Created on March 10, 2016, 12:07 AM
 *
 * This short intro added August 1, 2017 when I realized I hadn't signed this.
 */

#include "MathMethods.h"
#include <stdlib.h>
#include <math.h>
//#include "Stack.h"

/******************************************************************************
 *                                                                            *
 *                    CONSTANTS, NUMBERS, AND VARIABLES                       *
 *                                                                            *
 ******************************************************************************/

//unless otherwise defined, DELTA is a small number near zero
static double DELTA = EPSILON;

/*This function is used to set a static system variable called DELTA.
 *DELTA is used to take integrals, derivatives, and other values
 */
void setDelta(double d)
{
	DELTA = d;
}

/******************************************************************************
 *                                                                            *
 *                       MEMORY MANAGEMENT OF MATRICES                        *
 *                                                                            *
 ******************************************************************************/
/**
 * function: initMat()
 * just a shortcut for doing calloc() on a matrix.
 *
 * All memory is managed a pointer-pointers. This function allocates all needed memory.
 *
 * 
 * @param m the number of rows of the matrix mat
 * 
 * @param n the number of cols of the matrix mat
 * 
 * @param mat the matrix you wish to initialize
 */
void initMat(int rows, int cols, double*** matPtr)
{
	int k = 0;
	*matPtr = (double **) calloc(rows ,sizeof(double*));
	for (k = 0; k < rows; k++) {
		(*matPtr)[k] = (double *) calloc(cols, sizeof(double));
	}

}

/**
 * function: initMat()
 * just a shortcut for freeing a matrix.
 *
 * All memory is managed a pointer-pointers. This function deallocates all used memory.
 *
 * 
 * @param m the number of rows of the matrix mat
 * 
 * @param mat the matrix you wish to initialize
 */
void freeMat(int rows, double*** matPtr)
{
	for(int k = 0; k < rows; k++){
		free((*matPtr)[k]);
	}
	free(*(matPtr));
}
/** 
 * Function: copyArry()
 * Copy the values of array "from" into array "to".
 * Just a helper function for copying from one 1D array to another
 * 
 * @param n the number of entries
 * 
 * @param from the array with desired values
 * 
 * @param to the array to copy into
 *
 */
static void copyArray(int n, double *from, double *to)
{	
	int j = 0;
	for (j = 0; j < n; j++) {
		to[j] = from[j];
	}
}

/******************************************************************************
 *                                                                            *
 *                         LINEAR ALGEBRA OPERATIONS                          *
 *                                                                            *
 ******************************************************************************/

/**
 * Function: getRow()
 * Fetch the values of the r number row and store the first column
 * of the argument argument "store"
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
void getRow(int m, int n, double** mat, int r, double** store)
{
	int j=0;
	for (; j < n; j++) {
		store[j][0] = mat[r][j];
	}
	return;
}

/**
 * Function: getCol()
 * Fetch the values of the c number col and store the in the first column
 * of the argument argument "store"
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
void getCol(int m, int n, double** mat, int c, double** store)
{
	int i = 0;
	for (i = 0; i < m; i++) {
		store[i][0] = mat[i][c];
	}
	return;
}

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
double dotProduct(int m, double** v1, double** v2)
{
	double dot = 0;
	int i = 0;
	for (i = 0; i < m; i++) {
		dot += v1[i][0] * v2[i][0];
	}
	return dot;
}

/** 
 * Function: euclideanNorm()
 * Compute the euclidean norm (aka magnitude or spatial length) of a given 
 * real valued m-vector
 * 
 * @param m the length of the vector v
 * 
 * @param v the vector which will be evaluated
 */
double euclideanNorm(int m, double** v)
{
	double squareOfMagnitude = 0;
	int j=0;
	for( j = 0; j < m; j++)
	{
		squareOfMagnitude += v[j][0] * v[j][0];
	}

	return sqrt(squareOfMagnitude);
}

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
void swapRows(int n, double** mat, int r1, int r2)
{
	int j=0;
	double temp;
	for (j = 0; j < n; j++) {
		temp = mat[r1][j];
		mat[r1][j] = mat[r2][j];
		mat[r2][j] = temp;
	}
	return;
}

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
void scaleRow(int n, double** mat, int r, double scaleBy)
{
	int j=0;
	for (j = 0; j < n; j++) {
		mat[r][j] *= scaleBy;
	}
	return;
}

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
void elimRow(int n, double** mat, int r1, int r2, double elimBy)
{
	int j=0;
	for (j = 0; j < n; j++) {
	mat[r1][j] -= elimBy * mat[r2][j];
	}
	return;
}

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
void copyMatrix(int m, int n, double** from, double** to)
{
	int i=0;
	int j=0;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			to[i][j] = from[i][j];
		}
	}
}
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
void transpose(int m, int n, double** mat, double** matTranspose)
{
	//iterators
	int i = 0;
	int j = 0;

	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			matTranspose[j][i] = mat[i][j];
		}
	}
}

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
 *		  memory of a matrix. mat will be modified by this algorithm, so only 
 *		  give this function a matrix you want to change.
 * 
 * @return if the matrix was square, return the determinant.
 *		  if the matrix was not square, we return a zero, which is just a 
 *		  placeholder value so that we dont return garbage.
 */
double ref(int m, int n, double** mat)
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
	//this next int flag is for if we find a row of zero (i.e. not invertible)
	double det = 1;
	//check if our matrix is square.
	if (m != n) {
		det = 0;
	}
	while (i < m && j < n) {
		/* due to floating point rounding, we may end up with the difference
		 * between a number and another being nonzero, even when by our
		 * calculations by hand, they are equal. The constant
		 * DELTA is a buffer to make sure we still treat zeros
		 * properly. A note: fabs(x) is a function such that it takes a real
		 * number and returns the real number's magnitude. "fabs" is short for
		 * "floating point absolute value"
		 */
		if (fabs(mat[i][j]) < DELTA) {//if matrix's ij value is near zero
			for (x = i; x < m; x++) {//iterate through rows (vertically)
				if (fabs(mat[x][j]) > DELTA) {
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
				//due to order of computation, this loop is computed backwards.
				//We do not use the elimRow() function because the lower limit
				//of the loop is j not necessarily zero.
				//this saves a bit of computation time
				mat[x][y] -= mat[i][y] * mat[x][j];
			}
		}
		j++;
		i++; //all good, iterate forward!
	}
	//all done, return the determinant
	return det;
}

/** 
 * Function: rref() AKA Reduced Row Echelon Form
 * This function calls ref() once to make sure that it toying with a matrix
 * in row echelon form. After that it starts iterating from row (m - 1) all 
 * the way to the top, eliminating the values of the leading ones as it goes.
 * 
 * Keep in mind that if you write:
 * 
 *	  double mat[m][n] = {...some set of vectors...};
 *	  ref(m, n, mat);
 *	  double det = rref(mat);
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
 *		  return zero. We assume you understand that this zero has no
 *		  correlation to a determinant for the m not equal to n case.
 */
double rref(int m, int n, double** mat)
{
	double det = ref(m, n, mat);
	//these are our loop iterators, x and y:
	int x=0;
	int y=0;
	//these are the values of the row/column we are working on
	int i = m - 1;
	int j = 0;
	//first find the first nonzero row from the bottom up:
	//this is the thing that tells us we found a pivot column
	int foundPivotFlag = 0;
	while (i >= 0) {
		//first we find the first pivot position, from the bottom up
		for (x = i; x >= 0; x--) {
			for (y = 0; y < n; y++) {
				if (mat[x][y] > DELTA) {
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
double invert(int m, int n, double** mat)
{
	double det = 0;
	int i=0;
	int j=0;

	if (m != n) {
		return 0;
	}
	
	
	double** checkMat;
	initMat(m, n, &checkMat);
	//use this matrix to check what the determinant is before we try to find
	//and inverse.

	copyMatrix(m, n, mat, checkMat);

	det = ref(m, n, checkMat);
	freeMat(m, &checkMat);
	if (fabs(det) < DELTA) {
		return 0;
	}
	
	double** augMat;
	initMat(m, 2*n, &augMat);
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
	freeMat(m, &augMat);
	return 1/det;

}

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
void matrixMultiply(int m, int n, int p, double** left, double** right,
		double** product)
{

	int i = 0;
	int j = 0;
	double** v1;
	double** v2;
	initMat(n, 1, &v1);
	initMat(n, 1, &v2);
	for (i = 0; i < m; i++) {
		for (j = 0; j < p; j++) {
			getRow(m, n, left, i, v1);
			getCol(n, p, right, j, v2);
			product[i][j] = dotProduct(n, v1, v2);
		}
	}

	freeMat(n, &v1);
	freeMat(n, &v2);
}

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
void arrayAdd(int m, double *left, double *right, double *sum)
{
	int j;
	for(j = 0; j < m; j++){
		sum[j] = left[j] + right[j];
	}
}
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
void arrayScale(int m, double *arr, double scaleBy, double *scaled)
{
	int j;
	for(j = 0; j < m; j++){
		scaled[j] =  scaleBy*arr[j];
	}
}

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
void matrixAdd(int m, int n, double** left, double** right, 
				double** sum)
{
	int j;
	for(j = 0; j < n; j++){
		arrayAdd(m, left[j], right[j], sum[j]);
	}
}

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
double rectIntegral(double (*funct)(double), double a, double b)
{
	//value of a reimann sum
	double riemann = 0;
	//index to calulate an integral with
	double index = 0;
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
		riemann += funct(index) * DELTA;
		index += DELTA;
	}

	//multiply by the factor mentioned to get the finite integral approximation
	return riemann * backwards;
}

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
double simpIntegral(double (*funct)(double), double a, double b)
{
	//value of a reimann sum
	double simpson = 0;
	//index to calulate an integral with
	long int index = 0;
	//n is the number of iterations
	long int n = (long int)abs((long int)((a - b)/DELTA));
	//a flag and swap value for when you integrate from higher to lower
	double backwards = 1;
	//x_k and x_kp1 are the k^th and (k+1)^th points in the integration domain
	double x_k = 0;
	double x_kp1 = 0;
	//x_bar is the average of x_k and x_kp1 at the k^th iteration
	double x_bar = 0;

	//use negative of an integral property
	if (b < a) {
		backwards = b;
		b = a;
		a = backwards;
		backwards = -1;
	}

	//first iteration, x_k is left bound.
	x_k = a;
	x_kp1 = a + DELTA;

	//calulate the riemann sum
	while (index < n) {
		x_bar = (x_k + x_kp1)/2;
		simpson += (funct(x_k) +  4*funct(x_bar) + funct(x_kp1))*DELTA;
		x_k += DELTA;
		x_kp1 += DELTA;
		//printf("index: %d %f\n", index, simpson);
		index++;
	}

	//multiply through by the factors needed for the simpson rule.
	return (simpson * backwards)/6;
}

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
double simpleDerivative(double (*funct)(double), double a)
{
	return (funct(a + DELTA) - funct(a)) / DELTA;
}

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
void adamsBash3(int n, void (*vecField)(double*, double*),
				double** y_initial, double deltaT, 
				int iostep, long int numPnts,
				double *(*y_solve),//[(numPnts - (numPnts % iostep))/iostep], 
				double *time)
{
	//Initialize all memory needed to solve ODE:
	double *tempVec1 = calloc(n, sizeof(double));
	double *tempVec2 = calloc(n, sizeof(double));
	double *yjm2 = calloc(n, sizeof(double));
	double *yjm1 = calloc(n, sizeof(double));
	double *yj = calloc(n, sizeof(double));
	double *yjp1 = calloc(n, sizeof(double));
	long int i = 0;	
	int j = 0;
	int numSaved = -1;
	double tj = 0;
	time[0] = tj;
	copyArray(n, y_initial[0], yjm2);

	//These Several lines are equivalent to the MATLAB code:
	//yjm1 = y0 + DT*0.5*(vecField(y0+DT*vecField(y0,)) + vecField(y0));
	vecField(y_initial[0], tempVec1);
 	arrayScale(n, tempVec1, deltaT, tempVec1);
	arrayAdd(n, y_initial[0], tempVec1, tempVec1);
	vecField(tempVec1, tempVec2);
	vecField(y_initial[0], tempVec1);
	arrayAdd(n, tempVec1, tempVec2, tempVec1);
	arrayScale(n, tempVec1, 0.5*deltaT, tempVec1);
	arrayAdd(n, y_initial[0], tempVec1, yjm1);
	//Step Time Forward
	tj += deltaT;
	if(1 == iostep){
		time[1] = tj;
	}
	
	//These lines are equivalent to the MATLAB code: (fun is vecField)
	//yj = yjm1 + DT*0.5*(fun(yjm1+DT*fun(yjm1)) + fun(yjm1));
	vecField(yjm1, tempVec1);
	arrayScale(n, tempVec1, deltaT, tempVec1);
	arrayAdd(n, yjm1, tempVec1, tempVec1);
	vecField(tempVec1, tempVec2);
	vecField(yjm1, tempVec1);
	arrayAdd(n, tempVec1, tempVec2, tempVec1);
	arrayScale(n, tempVec1, 0.5*deltaT, tempVec1);
	arrayAdd(n, tempVec1, yjm1, yj);
	//Step Time Forward
	tj += deltaT;
	//if( 2 % iostep == 0){
	// 	time[2] = tj;
	//}
	//Store the current solution values:
	for (j = 0; j < n; j++) {
		y_solve[j][0] = yjm2[j];
		numSaved++;
	}
	if(1 == iostep){
		for (j = 0; j < n; j++) {
			y_solve[j][numSaved] = yjm1[j];
		}
		numSaved++;
	}
	if(2 % iostep == 0){
		for (j = 0; j < n; j++) {
			y_solve[j][numSaved] = yj[j];
		}
		time[numSaved] = tj;
		numSaved++;
	}

	//Done with kickstart, starting the AB3 iterations
	i = 3;
	while(i < numPnts){
		//Step Time Forward
		tj += deltaT;
		// The following section is equivalent to the MATLAB code:
		//yjp1= yj+(DT/12)*(23*fun(yj)-16*fun(yjm1)+5*fun(yjm2));
		vecField(yj, tempVec1);
		arrayScale(n, tempVec1, 23, tempVec1);
		vecField(yjm1, tempVec2);
		arrayScale(n, tempVec2, -16, tempVec2);
		arrayAdd(n, tempVec1, tempVec2, tempVec1);
		vecField(yjm2, tempVec2);
		arrayScale(n, tempVec2, 5, tempVec2);
		arrayAdd(n, tempVec1, tempVec2, tempVec1);
		arrayScale(n, tempVec1, deltaT/12.0, tempVec1);
		arrayAdd(n, yj, tempVec1, yjp1);
		//Every iostep iterations, we must save the state of the ODE at time tj:
		if( (i % iostep == 0) &&
				(numSaved < (numPnts - (numPnts % iostep))/iostep)){
			for (j = 0; j < n; j++) {
				y_solve[j][numSaved] = yjp1[j];
			}
			time[numSaved] = tj;
			numSaved++;
		}
		//Iteration i done. Now shifting old values over by 1 step.
		copyArray(n, yjm1, yjm2);
		copyArray(n, yj, yjm1);
		copyArray(n, yjp1, yj);
		i++;
	}
	//This last section does some handling for when someone sets iostep=numPnts
	//Basically, we choose to give back the last iteration of the solve.
	//numPnts == 1 means only return inital condition.
	if(iostep == numPnts){
		if(iostep == 1){//numPnts == 1 means only return inital condition.
			for (j = 0; j < n; j++) {
					y_solve[j][0] = yjm2[j];
				}
			time[0] = 0;
		}
		if(iostep == 2){//numPnts == 2 means only return 1 step forward
			for (j = 0; j < n; j++) {
					y_solve[j][0] = yjm1[j];
				}
			time[0] = deltaT;
		}
		if(iostep == 3){//numPnts == 3 means only return 2 steps foward.
			for (j = 0; j < n; j++) {
					y_solve[j][0] = yj[j];
				}
			time[0] = deltaT + deltaT;
		}
		if(iostep > 3){//numPnts == 1 means only return late position
			for (j = 0; j < n; j++) {
					y_solve[j][0] = yjp1[j];
				}
			time[0] = tj;
		}
	}
	free(tempVec1);
	free(tempVec2);
	free(yj);
	free(yjm1);
	free(yjm2);
	free(yjp1);
}

/**
 * Function: binomialCoef()
 * Finds the "n choose r" value. also known as the binomial coefficient.
 *
 * @param n the size of the set to choose from
 * @param r the number choose from a set.
 */
long int binomialCoef( long int n, long int r )
{
	//loop iterator
	int i = 0;
	//the return value
	int pascal = n;
	if (r > n) {
		return 0; //invalid choice case
	}

	if ((r * 2) > n) {
		//the binomial coefficient is symmetric about the center of
		//Pascal's triangle, so we use this one weird trick:
		r = n-r; //click here to find out my one trick! mathematicians hate me!
	}

	if (r == 0) {
		return 1; //only 1 way to choose nothing
	}

	for(i = 2; i <= r; i++ ) { //this simulates traversing Pascal's triangle
		pascal *= (n-i+1);
		pascal /= i;
	}
	//it is only fitting that pascal would traverse his triangle.
	return pascal;
}

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
void vanderInterpolate(long int n, double* x, double* y, double* c)
{
	//loop iterators

	int i = 0;
	int j = 0;

	//n is the number of coefficients. it is also the length of c, x, and y
	//this is essentially the most important number to the algorithm

	//matrix of x values raised to powers
	double** xMat;
	initMat(n, n, &xMat);
	//intialize xMat right after declaring
	//(to avoid eldritch C memory issues.)
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			xMat[i][j] = pow(x[i] , j);
			//value of of the (i,j) entry is the i number to the j power
		}
	}

	//column vector of y values
	double** yMat;
	initMat(n, 1, &yMat);
	//initialize yMat right after declaring
	//(to avoid eldritch C memory issues.)
	//This may seem odd, but it is done so we can use the linear algebra
	//functions written above. y[j] is technically a row vector.
	//we need a column vector to do the operations.
	for(i = 0; i < n; i++) {
		yMat[i][0] = y[i];
	}

	//column vector which is the product of inverse(xMat) and yMat
	double** pMat;
	initMat(n, 1, &pMat);
	//intialize pMat right after declaring
	//(to avoid eldritch C memory issues.)
	for(i = 0; i < n; i++) {
		pMat[i][0] = 0;
	}

	//now building a single polynomial:
	//take the inverse of the xMat matrix so we can find the
	//coefficients of our polynomial

	invert(n, n, xMat);
	//multiply this inverse matrix by yMat to get  one of the
	//polynomial coefficient vectors

	matrixMultiply(n, n, 1, xMat, yMat, pMat);

	//now that we have the coefficients of our polynomial,
	//load the output into c.
	for(i = 0; i < n; i++) {
		c[i] = pMat[i][0];
	}
	freeMat(n, &xMat);
	freeMat(n, &yMat);
	freeMat(n, &pMat);
	//done
	return;
}


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
/*Not Yet implemented
void lagrangeInterpolate(long int n, double x[n], double y[n], double c[n])
{
}*/

/**
 * Function: meanPolynomial()
 * This function will make a "best-fit" polynomial of a specified degree. Given
 * 2 vectors representing the x and y coordinates of a graph, the function will
 * analyze them in the following manner:
 *
 *		1) ->  Call the interpolate function using the given x and y vectors
 *
 *		2) -> Do this for recursively for every n size subset of our data 
 *			  points, taking the mean of the results as you compute. The 
 *			  result is a vector of coefficients for a best fit polynomial of 
 *			  degree d.
 *  
 *  (NOTE)  2)-> Recursive computation finds a num_coefs size subset of the x,y
 *			   data pairs then calls vanderInterpolate() on them to compute.
 *			   only guaranteed to work if you pass -1 top path 
 *				and pass NULL to *combo.
 *
 * @param path an iterator for recursive computation. 
 *				ALWAYS PASS AS -1 TO INITIALIZE
 * 
 * @param *combo this is for combinatorially iterating through possible 
 *				  indeces when interpolating. 
 *				  ALWAYS PASS AS NULL POINTER TO INITIALIZE!
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
 *			   the polynomial. c[0] is the constant term, c[d] is the highest
 *			   power's coefficient.
 */
void meanPolynomial(long int path, Stack *combo, long int num_points, 
					double x[num_points], double y[num_points], 
					long int num_coefs, double c[num_coefs])
{   	
	//loop iterator
	long int i = 0;
	//these two arrays hold the particular data for a given subset of the points
	double combo_x[num_coefs];
	double combo_y[num_coefs];
	//this array is similar to the two above, but for coefficients
	double combo_c[num_coefs];
	long int popValue = 0;
	if(path == -1) {
		//Initialize Stack at the start of the recursion
		StackInit(&combo, num_coefs);
	} else {
		//already partially recursed
		StackPush(combo, path);
	}
	//if we have a subset (not end of path)
	if(StackGetSize(combo) == num_coefs) {
		//process combinatorial data
		for(i = 0; i < num_coefs; i++) {
			//make the placeholder arrays equal to the specific x, y values
			combo_x[i] = x[combo->stackItems[i]];
			combo_y[i] = y[combo->stackItems[i]];
		}
		vanderInterpolate(num_coefs, combo_x, combo_y, combo_c);
		for(i = 0; i < num_coefs; i++) {
			c[i] += combo_c[i];
		}
		StackPop(combo, &popValue);
		return;
	}
	for (i = path + 1; i < num_points; i++) { //go down all paths > current
		meanPolynomial(i, combo, num_points, x, y, num_coefs, c);
	}
	StackPop(combo, &popValue);
	if(path == -1) {
		StackFree(combo);

		for(i = 0; i < num_coefs; i++) {
			c[i] /= binomialCoef(num_points, num_coefs);
		}
	}
	return;
}

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
 *
 * @param n the number of input data points
 *
 * @param x[n] array of x-axis data points
 *
 * @param y[n] array of y-axis data points
 *
 * @param m the number of functions
 *
 * @param (*phi[m])(double) array of pointers to functions which take a
 * 							double and return a double.
 *
 * @param result[m] array of coeficients to functions (*phi[m])(double)
 */
void discreteLeastSquares(int n, double* x, double* y, int m,
							double (*(*phi))(double), double* result)
{
	int i = 0;
	int j = 0;
	//matrix of phi functions evaluated at x values
	double *(*bMat);
	initMat(n, m, &bMat);
	//the transpose of bMat
	double *(*bMatTr);
	initMat(m, n, &bMatTr);
	/* bMat is n by m, which goes against the usual convention of alphabetical
	*  order of the size. I do this though because I prefer the convention
	*  of letting n represent the data set, as in previous function fit models.
	*/
	//This matrix is the n by n+1 which will have bMatTr*bMat in first n by n
	//entries. In the last column, we have bMatTr*y. We do rref() on this matrix
	//and the result in the last column is the least squares coordinates.
	double *(*leastSqr);
	initMat(m, m+1, &leastSqr);
	//This n by n matrix has the sole purpose of being the product of
	//and bMatTr. It will be copied into leastSqr.
	double *(*tempMat);
	initMat(m, m, &tempMat);
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			bMat[i][j] = (*phi[j])(x[i]);
		}
	}
	transpose(n, m, bMat, bMatTr);
	matrixMultiply(m, n, m, bMatTr, bMat, tempMat);
	//Now multiply bMatTr and y. Not using matrix multiply because y is
	//in the wrong format. We also load tempMat into leastSqr.
	for(i = 0; i < m; i++){
		for(j = 0; j < n; j++){
			leastSqr[i][m] += bMatTr[i][j] * y[j];
			if(j < m){
				leastSqr[i][j] = tempMat[i][j];
			}
		}
	}
	//solve the system of equations arrising from the least squares condition
	rref(m, m+1, leastSqr);
	for(i=0; i < m; i++){
		result[i] = leastSqr[i][m];
	}
	freeMat(n, &bMat);
	freeMat(m, &bMatTr);
	freeMat(m, &leastSqr);
	freeMat(m, &tempMat);
}

double discreteFourier(int n, double x[n], double y[n], int numFreqs, 
						double cosines[numFreqs], double sines[numFreqs])
{
	int i = 0;
	int j = 0;
	register int TwiceNumFreqs = 2*numFreqs;
	register double angularFreq = 2*M_PI/(x[n-1] - x[0]);
	//This double contains the constant term of the fourier transform
	double meanOfData = 0;
	//matrix of phi functions evaluated at x values
	double *(*bMat);
	initMat(n, TwiceNumFreqs, &bMat);
	//the transpose of bMat
	double *(*bMatTr);
	initMat(TwiceNumFreqs, n, &bMatTr);
	//This matrix is the n by n+1 which will have bMatTr*bMat in first n by n
	//entries. In the last column, we have bMatTr*y. We do rref() on this matrix
	//and the result in the last column is the least squares coordinates.
	double *(*leastSqr);
	initMat(TwiceNumFreqs, (TwiceNumFreqs)+1, &leastSqr);
	//This n by n matrix has the sole purpose of being the product of
	//and bMatTr. It will be copied into leastSqr.
	double *(*tempMat);
	initMat(TwiceNumFreqs, TwiceNumFreqs, &tempMat);
	for(i = 0; i < n; i++){
		meanOfData += y[i];
		for(j = 0; j < numFreqs; j++){
			bMat[i][2*j] = cos(angularFreq * (j+1) * (x[i]));
			bMat[i][(2*j)+1] = sin(angularFreq *(j+1) * (x[i]));
		}
	}
	transpose(n, TwiceNumFreqs, bMat, bMatTr);
	matrixMultiply(TwiceNumFreqs, n, TwiceNumFreqs, bMatTr, bMat, tempMat);
	//Now multiply bMatTr and y. Not using matrix multiply because y is
	//in the wrong format. We also load tempMat into leastSqr.
	for(i = 0; i < TwiceNumFreqs; i++){
		for(j = 0; j < n; j++){
			leastSqr[i][TwiceNumFreqs] += bMatTr[i][j] * y[j];
			if(j < TwiceNumFreqs){
				leastSqr[i][j] = tempMat[i][j];
			}
		}
	}
	//solve the system of equations arrising from the least squares condition
	rref((TwiceNumFreqs), (TwiceNumFreqs)+1, leastSqr);
	for(i=0; i < numFreqs; i++){
		cosines[i] = leastSqr[2*i][TwiceNumFreqs];
		sines[i] = leastSqr[(2*i)+1][TwiceNumFreqs];
	}
	freeMat(n, &bMat);
	freeMat(TwiceNumFreqs, &bMatTr);
	freeMat(TwiceNumFreqs, &leastSqr);
	freeMat(TwiceNumFreqs, &tempMat);
	return meanOfData / n;
}
