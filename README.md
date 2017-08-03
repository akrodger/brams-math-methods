# brams-math-methods
A collection of math computation C libraries written by Abram "Bram" Rodgers.
README by Abram Rodgers
# Contact: AKS.Rodgers (at) gmail (dot) com
-------------
| Contents: |
-------------
- Intro: Origin of this Library
- Section 0: Library Usage
- Section 1: Linear Algebra
- Section 2: Calculus and Numerical Methods
- Section 4: Integer Stack Implementation
- Section 5: Special Thanks

---------------------------------
| Intro: Origin of this Library |
---------------------------------
This program was initially the primary content for a school report. When I was a student in MATH 181 History of Mathematics (Winter 2016) at UC Santa Cruz, we all had to do a historical report on a famous mathematician. We also had to showcase some of that mathematician's work. I chose to write about Carl Gauss and showcase Gaussian Elimination, a particular matrix algorithm which can find a solution for systems of linear equations of nth order. Since this is a pretty standard exercise for students in a computer science course, I figured I could probably write my own version of the very famous algorithm for my math course. After approval from the professor, I wrote up the code over the course of a week or so and included it in my historical report with gratuitous comments and documentation (since he did not know C). The school project was very successful and writing more in the library became more of a side project.

Today the library sits on GitHub to show some of my personal coding skills. Since no algorithm in this package was particularly novel, simply my own take on various computational problems, I decided to make it all public domain at first. Plenty of alternate versions of the Linear Algebra algorithms are open source. However, after writing the numerical methods functions, I decided to open source the whole package, since statistical software is a bit more novel than some intro Linear Algebra algorithms.

If you wish to download and use my code in particular, I would be curious to know what it is being used for. If you wish to ask a question, collaborate on some project, or something else, then just shoot me an email.


----------------------------
| Section 0: Library Usage |
----------------------------
Program located in folder bmm.
Library Dependencies: 
	- gcc
	- make
	- gnuplot

Usage:

- Run make in the bmm directory

- LinAlg : This program is a showcase of the Linear Algebra operations. It loads a matrix with values that are known to be invertible (the i-j entry is just i to the power of j). It then inverts this matrix and multiplies the inverse by the original matrix to confirm that a matrix times its inverse is the identity. run ./LinAlg for usage instructions. (You get to choose the dimensions. Note that if you give it a non-square matrix, the inverse has no actual meaning, nor does the resultant multiplication)

- interpolator : This program reads in a list of (x,y) data pairs from the file "bmm/fio/InterpolatorPoints.txt". Then it will run meanInterpolate() on those points.

InterpolatorPoints.txt formatting:

    - Any line starting with a '~' is skipped. Any text after that is not read at all.

    - A line starting with a "t:" (colon optional) indicates that the line has an integer which indicates how many terms are in the polynomial you want to approximate. (Terms in a polynomial is equal to the degree plus 1.)

    - A line starting with 'x' indicates that the line contains the x coordinates of the array of data points serparated by spaces. There may not be trailing sapces after the last number..

    - A line starting with 'y' indicates that the line contains the y coordinates of the array of data points serparated by spaces. There may not be trailing sapces after the last number.

    - A line starthing with any other character throws and error and tells you to format the file properly.

    - Parser also ignores a line which contains no text. That is, it is just nothing but a "newline" or "line feed". This is treated the same as though it saw a '~'.

    - Note: No lines may between 'x' or 'y' and their corresponding coordinates on the next line.

Example of a valid file:

    ~ First line of file
    ~ "This text is not read by the parser"
    t: 5
    x: The colon is not read by the parser, only the x. This text also not read.
    1   2   4   6    7
    ~ As many spaces as you want can be between numbers

    ~ The above line is ignored completely.
    y:
    3 4 -6 6     -8
    ~ Last Line of file

- OdeSovler : This program uses the Adams-Bashforth ODE solver to compute a trajector of an LpCircle. Plotting is done in gnuplot using freely available code by N. Devillard. (Devillard's code included in bmm/src folder.)

-----------------------------
| Section 1: Linear Algebra |
-----------------------------
This library contains several algorithms from an introductory Linear Algebra course.

Notable functions are:
- dot product
- row swap
- row scale
- row elimination
- copy matrix
- transpose matrix
- row echelon form (Gaussian Elimination)
- reduced row echelon form
- invert matrix
- multiply matrices

Planned Functions:
- QR Factorization.
Currently, the design problem with QR is how to represent complex eigenvalues. I'm considering forking this library to C++ and using a complex numbers class. Considering that using complex numbers would minimally double the number of floating point operations per function. I'm not sure if this is a good idea. I may work on it just for the mathematical exercise though.

---------------------------------------------
| Section 2: Calculus and Numerical Methods |
---------------------------------------------
There are several functions here:
- Two types of finite integral with small delta approximation
- Derivative with small delta approximation
- Adams-Bashforth Orinary 3 step Differential Equation Solver
- Vandermonde polynomial interpolation
- Mean Polynomial approximation
- Least Squares Approximation (Unstable, still in testing)

Adams-Bashforth is a method for using local polynomial interpolation to find the solutions to n-dimensional ODEs. I implement it here in full generality. Keep note that this function may use a ton of memory if you set the IO step too low.

Mean polynomial approximation takes a set of (x, y) data points and forms a polynomial with them of a specified degree. It makes use of the linear algebra operations defined in the library to do this. Read the header file for more information.

Least Squares is an algorithm commonly used in statistical analysis to find a function which represents a set of measurements. We give the algorithm a linearly independent array of pointers to functions of the form 

	double foo(double)

We then compute the projection matrix of the x data points with respect to this array, minimizing the error from a linear combination of functions and the y values. The result is a set of coordinates, equivalent to the coeficients of either of polynomial approximation functions mentioned above.


-------------------------------------------
| Section 3: Integer Stack Implementation |
-------------------------------------------
This is an implementation of a stack data structure written by me. The reason for this is simply because having a stack to run computations with is useful. Currently, it is only used in the Numerical Methods section. MathMethods depends on Stack, so make sure to have them in the same directory.

-----------------------------
| Section 4: Special Thanks |
-----------------------------
- James Iwamasa

In fall 2016, James Iwamasa, a personal friend from the computer science department at UC Santa Cruz, helped me implement the meanPolynomial() function. For the mean polynomial approximation  algorithm, James implemented a recursive algorithm in C++ using vector objects to find all subsets of a list of integers. Then I converted this code to C, restructured it to work with the stack library, then optimized memory usage so that my laptop wouldn't crash due to RAM filling up for large enough sets of integers.
