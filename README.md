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

Today the library sits on GitHub to show some of my personal coding skills. Since no algorithm in this package is particularly novel, simply my own take on various computational problems, I decided to make it all public domain. Plenty of alternate versions of these algorithms are open source. If you wish to download and use my code in particular, I would be curious to know what it is being used for. If you wish to ask a question, collaborate on some project, or something else, then just shoot me an email.

Run make in brams-math-methods to compile program. Note, you make run into an error noting that a file named 'obj' does not exist. If you enounter this, make obj in the same directory as the makefile and try to make again.

----------------------------
| Section 0: Library Usage |
----------------------------
Library Dependencies: Nothing other than gcc and make

Usage:

- Run make in the bmm directory

- LinAlg : This program is a hard coded demonstration of the linear Algebra operations

- Numer : This is a nice little Interpolation Tester Program. Prints out a polynomial of the form (n * x^(n-1)) (for the nth term) after calculating a number of points to interpolate along. Uses the meanInterpolate() function, a combinatorial version of Lagrange Interpolation. Usage instructions built into the program it, run ./Numer for instructions.

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
