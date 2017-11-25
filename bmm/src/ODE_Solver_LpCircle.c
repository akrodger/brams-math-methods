// **** Include libraries here ****
// Standard libraries
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
// User libraries
#include "MathMethods.h"
#include "gnuplot_i.h"

#ifndef NULL
#define NULL 0
#endif

#define ERR_MSG \
	"\nUsage: OdeSolve [P] [Num Pnts] [IO Step] [Delta T] [A] [B] [W]\n\n"\
	"P:         An even integer. This is the P in an LP norm used in ODE.\n"\
	"Num Pnts:  A long int. This is the number of iterations for AB3.\n"\
	"IO Step:  An int. This is number of iterations between saving.\n"\
	"Delta T:   A double. This is the time step of the ODE solver.\n"\
	"A:         A double. This is x-direction maximum of Steady State.\n"	\
	"B:         A double. This is y-direction maximum of Steady State.\n"	\
	"W:         A double. This is Time Scaling of Steady State.\n"\
	"                     For A=B=1, it is the angular frequency.\n\n"\
	"Default: A = B =  W = 1, Num Pnts = 1000000, int IO Step = 1, and " \
	"Delta T =  8e-5. All arguments other than P are optional. P must be even."\
	"\n\n"\
// **** Set macros and preprocessor directives ****


// **** Define global, module-level, or external variables here ****

int P;
double A;
double B;
double W;
// **** Declare function prototypes ****
void LP_OSC(double arg[2], double ret[2])
{
	ret[0] = W*((arg[0]/A)*(1 - (pow(arg[0]/A,P)  +
				 pow(arg[1]/B,P))) - pow(arg[1]/B,P-1));
	ret[1] = W*((arg[1]/B)*(1 - (pow(arg[0]/A,P)  +
				 pow(arg[1]/B,P))) + pow(arg[0]/A,P-1));
}


int main(int argc, char **argv)
{	
	
	int j;
	P = 2;
	long int NUM_PNTS = 1000000;
	int IO_STEP = 1;
	double DT =  8e-5;
	A  = 1;
	B  = 1;
	W  = 1;
	 if(argc == 1 || argc > 8){
		printf(ERR_MSG);
		return 1;
	}

	char *end = NULL;
	P = strtol(argv[1], &end, 0);
	if(P % 2 == 1){
		printf(ERR_MSG);
		return 1;
	}

	if(argc >= 3){
		NUM_PNTS = strtol(argv[2], &end, 0);
		if(argc >= 4){
			IO_STEP = strtol(argv[3], &end, 0);
			if(argc >= 5){
				DT = strtod(argv[4], &end);
				if(argc >= 6){
					A = strtod(argv[5], &end);
					if(argc >= 7){
						B = strtod(argv[6], &end);
						if(argc >= 8)
							W = strtod(argv[7], &end);
					}
				}
			}
		}
	}
	
	
	double** y0;
	initMat(2, 1, &y0);
	
	y0[0][0] = 1;
	y0[1][0] = 0;

	gnuplot_ctrl *h1 = gnuplot_init() ;
	
	long int memLevel = (NUM_PNTS - (NUM_PNTS % IO_STEP))/IO_STEP;
	
	double *t = (double*) calloc(memLevel, sizeof(double));
	double** yPlot;
	initMat(2, memLevel, &yPlot);
	
	double *cosineWave = (double*) calloc(memLevel, sizeof(double));
	double *sineWave = (double*) calloc(memLevel, sizeof(double));
	
	printf("\nBeginning Computation of ODE...\n");
	adamsBash3(2, &LP_OSC, y0, DT, IO_STEP, NUM_PNTS, yPlot, t);
	
	printf("\nDone computing ODE. Writing plots now...\n");
	for(j = 0; j < memLevel; j++)
	{
		cosineWave[j] = yPlot[0][j];
		sineWave[j] = yPlot[1][j];
		//printf("\n\ncosine: %f\nsine: %f\ntime: %f\n", 
		//			yPlot[0][j],yPlot[1][j], t[j]);
	}
    gnuplot_resetplot(h1) ;
    gnuplot_setstyle(h1, (char*)"linespoints") ;
	gnuplot_cmd(h1, "set terminal png");
	gnuplot_cmd(h1, "set output \"Lp_cosine.png\"");
    gnuplot_plot_xy(h1, t, cosineWave, memLevel, (char*)"Lp Cosine");
	sleep(2);
	printf("\nDone with cosine. Saved to \"Lp_cosine.png\"\n");
	gnuplot_cmd(h1, "set output \"Lp_sine.png\"");
    gnuplot_resetplot(h1) ;
    gnuplot_plot_xy(h1, t, sineWave, memLevel, (char*) "Lp Sine");
	sleep(2);
	printf("\nDone with sine. Saved to \"Lp_sine.png\"\n");
	
	
    gnuplot_resetplot(h1) ;
	gnuplot_cmd(h1, "set output \"Lp_Circ.png\"");
    gnuplot_plot_xy(h1, cosineWave, sineWave, memLevel, (char*)"Lp Circle");
	
	
	gnuplot_close(h1) ;
	printf("\nDone with circle. Saved to \"Lp_Circ.png\"\n\n");
	
	freeMat(2, &y0);
	freeMat(2, &yPlot);
	free(t);
	free(cosineWave);
	free(sineWave);
}
