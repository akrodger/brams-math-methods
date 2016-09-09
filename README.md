# brams-math-models
A collection of math computation C libraries written by Abram "Bram" Rodgers.
README by Abram Rodgers
Contact: AKS.Rodgers (at) gmail (dot) com
-------------
| Contents: |
-------------
Intro Origen of this Library
Section 1: Linear Algebra
Section 2: Calculus

Intro: Origen of this Library
This program was initially the primary content for a school report. When I was a student in MATH 181 History of Mathematics (Winter 2016) at UC Santa Cruz, we all had to do a historical report on a famous mathematician. We also had to showcase some of that mathematician's work. I chose to write about Carl Gauss and showcase Gaussian Elimination, a particular matrix algorithm which can find a solution for systems of linear equations of nth order. Since this is a pretty standard exercise for students in a computer science course, I figured I could probably write my own version of the very famous algorithm for my math course. After approval from the professor, I wrote up the code over the course of a week or so and included it in my historical report with gratutious comments and documentation (since he did not know C). The school project was very successful and writing more in the library became more of a side project.

Today the library sits on GitHub to show some of personal coding skills. Since no algorithm in this package is particularly novel, simply my own take on various computational problems, I decided to make it all public domain. Plenty of alternate versions of these algorithms are open source. If you wish to download and use my code in particular, I would be curious to know what it is being used for. If you wish to ask a question, collaborate on some project, or something else, then just shoot me an email.

Editor: Netbeans For Linux

-----------------------------
| Section 1: Linear Algebra |
-----------------------------
This library contains several algorithms from an introductory Linear Algebra course.

Notable fucntions are:
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

More is planned to be written once I have completed a QR factorization implementation.
