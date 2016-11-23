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
- Section 2: Calculus
- Section 3: Numerical Methods
- Section 4: Integer Stack Implementation

---------------------------------
| Intro: Origin of this Library |
---------------------------------
This program was initially the primary content for a school report. When I was a student in MATH 181 History of Mathematics (Winter 2016) at UC Santa Cruz, we all had to do a historical report on a famous mathematician. We also had to showcase some of that mathematician's work. I chose to write about Carl Gauss and showcase Gaussian Elimination, a particular matrix algorithm which can find a solution for systems of linear equations of nth order. Since this is a pretty standard exercise for students in a computer science course, I figured I could probably write my own version of the very famous algorithm for my math course. After approval from the professor, I wrote up the code over the course of a week or so and included it in my historical report with gratuitous comments and documentation (since he did not know C). The school project was very successful and writing more in the library became more of a side project.

Today the library sits on GitHub to show some of my personal coding skills. Since no algorithm in this package was particularly novel, simply my own take on various computational problems, I decided to make it all public domain at first. Plenty of alternate versions of the Linear Algebra algorithms are open source. However, after writing the numerical methods functions, I decided to open source the whole package, since statistical software is a bit more novel than some intro Linear Algebra algorithms.

If you wish to download and use my code in particular, I would be curious to know what it is being used for. If you wish to ask a question, collaborate on some project, or something else, then just shoot me an email.


----------------------------
| Section 0: Library Usage |
----------------------------
Library Dependencies: Nothing other than gcc and make. Program located in folder bmm.

Usage:

- Run make in the bmm directory

- LinAlg : This program is a showcase of the Linear Algebra operations. It loads a matrix with values that are known to be invertible (the i-j entry is just i to the power of j). It then inverts this matrix and multipliea the inverse by the original matrix to confirm that a matrix times its inverse is the identity. run ./LinAlg for usage instructions. (You get to choose the dimensions. Note that if you give it a non-square matrix, the inverse has no actual meaning, nor does the resultant multiplication)

- interpolator : This program reads in a list of (x,y) data pairs from the file "bmm/fio/InterpolatorPoints.txt". Then it will run meanInterpolate() on those points.

InterpolatorPoints.txt formatting:

    - Any line starting with a '~' is skipped. Any text after that is not read at all.

    - A line starting with a "t:" (colon optional) indicates that the line has an integer which indicates how many terms are in the polynomial you want to approximate. (Terms in a polynomial is equal to the degree minus 1.)

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
- row echelon form (Gaussian Elimination)
- reduced row echelon form
- invert matrix
- multiply matrices

Planned Functions:
- QR Factorization

-----------------------
| Section 2: Calculus |
-----------------------
So far this only has 2 functions.
- finite integral with small delta approximation
- derivative with small delta approximation


--------------------------------
| Section 3: Numerical Methods |
--------------------------------
Two major operations are here:
- polynomial interpolation
- mean polynomial interpolation

This function takes a set of (x, y) data points and forms a polynomial with them of a specified degree. It makes use of the linear algebra operations defined in the library to do this. Read the header file for more information.

The intrpolate() and meanInterpolate() functions are considered UNSTABLE! Seems not to work very well for polynomials of degree hgiher than about 10. This could be due to a number of factors, including something possibly wrong in the linear algebra section.


More is planned to be written once I have completed a QR factorization implementation.
NOTE: QR factorization is shelved for now in favor of polynomial interpolation.

-------------------------------------------
| Section 4: Integer Stack Implementation |
-------------------------------------------
This is an implementation of a stack data structure written by me. The reason for this is simply because having a stack to run computations with is insanely useful. Currently, it is only used in the Numerical Methods section. MathMethods depends on Stack, so make sure to have them in the same directory.
