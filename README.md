# brams-math-methods
A collection of math computation C libraries written by Abram "Bram" Rodgers.
README by Abram Rodgers
# Contact: AKS.Rodgers (at) gmail (dot) com
-------------
| Contents: |
-------------
- Intro: What is this?
- Section 0: Least Square Regressional Analysis
- Section 1: Fourier Analysis Via Least Squares

------------------------
| Intro: What is this? |
------------------------
This is a branch of the open source mathematical modeling package titled "bmm". The goal of this branch is to show that the algorithms developed and maintained on the "master" branch are lightweight enough to run on an old microcontroller. One of the most well documented, well similulated, and well known ARM microcontroller systems is the Nintendo GameBoy Advance. Here, we use that system to demonstrate visually what the "bmm" library is capable of.

- I would like to highlight this note: None of the libraries used are dependent on the gameboy. You can use the code provided on the "master" branch to run on any system with a modern C compiler. (Any modern gcc should be fine.)

Inside this branch, you will find 2 tech demos which can be run on a Nintendo GameBoy Advance (referred to as GBA from here on) or on a GBA emulator.

Q: How do I run it?
	
	A:
	1) Download a gameboy advance emulator for your computer. I recommend mGBA ( found at https://mgba.io/ ) since it is cross platform, easy to set up, and open source. (Note: personally I use mednafen on Linux. Doesn't really matter which you use.)
	2) Go to the "LeastSquare" or "FourierTransform" folder.
	3) Run the ".gba" file using your emulator.

Q: How do I compile it?

	A:
	1) Download devkitPro, a GBA development kit ( found at https://devkitpro.org/ ).
	2) Make sure you have make installed (You probably installed it for devkitPro, my linux distro comes with it.)
	3) run "make clean" then "make" in either of the "LeastSquare" or "FourierTransform" folders.

-------------------------------------------------
| Section 0: Least Square Regressional Analysis |
-------------------------------------------------
This program generates a scatterplot using a trigonometric series with white noise added in. Then, it uses a 15 term polynomial to approximate that trig series. (See video.) The method used is the well known finite dimensional least squares regressive formula. (See derivation.)

-------------------------------------------------
| Section 1: Fourier Analysis Via Least Squares |
-------------------------------------------------
By setting the functions used in Least Squares to be trig functions, we actually can compute Fourier Coeficients of a scatterplot. The algorithm is the same, however we dynamically compute the trig functions used instead of requiring a set of function pointers. For details, see the header file in the main branch. (For visual, see video.)

