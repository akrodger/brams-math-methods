// **** Include libraries here ****
// Standard libraries
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// User libraries
#include "MathMethods.h"

#ifndef NULL
#define NULL 0
#endif
// **** Set macros and preprocessor directives ****

#define  ROWS 4
#define  COLS 4

// **** Define global, module-level, or external variables here ****

// **** Declare function prototypes ****

int main(int argc, char **argv)
{
	if(argc < 3 || 4 < argc)
	{
		printf("\nUsage: Numer [NUM_POINTS] [TERMS] [BUFFER]\n"
		"NOTE: 3rd argument is optional floating point percision setting.\n\n");
		return 1;
	}
	
	char *end = NULL;
	long int NUM_POINTS = strtol(argv[1], &end, 0);
	long int TERMS = strtol(argv[2], &end, 0);
	
	if(argv[3] != NULL)
	{
		double deltaBuffer = strtod(argv[3], &end);
		setDelta(deltaBuffer);
	}
	
    int j = 0;
	int i = 0;

    double *x = calloc(NUM_POINTS, sizeof(double));
    
    for(j = 0; j < NUM_POINTS; j++)
    {
        x[j] = j;
    }
    /*   
    x[0]  = 276;
    x[1]  = 278;
    x[2]  = 280;
    x[3]  = 282;
    x[4]  = 284;
    x[5]  = 286;
    x[6]  = 288;
    x[7]  = 290;
    x[8]  = 292;
    x[9]  = 294;    
    x[10] = 296;   
    x[11] = 298;   
    x[12] = 300;   
    x[13] = 302;   
    x[14] = 304;   
    x[15] = 306;   
    x[16] = 308;   
    x[17] = 310;  
    x[18] = 312;
    */
    double *y = calloc(NUM_POINTS, sizeof(double));
    
    
    for(j = 0; j < NUM_POINTS; j++)
    {
		for(i = 0; i < TERMS; i++)        
		{		
			y[j] += (i+1) * pow((x[j]), i);
		}
    }
     /* 
    y[0]  = 2.65;    
    y[1]  = 27;
    y[2]  = 42;
    y[3]  = 55.3;
    y[4]  = 63;
    y[5]  = 73;
    y[6]  = 78;
    y[7]  = 81;
    y[8]  = 84;
    y[9]  = 85;    
    y[10] = 83;   
    y[11] = 78;   
    y[12] = 74;   
    y[13] = 69;   
    y[14] = 62;   
    y[15] = 56;   
    y[16] = 39;   
    y[17] = 16;  
    y[18] = 1.27;
    */
    
	
    double *arr = calloc(NUM_POINTS, sizeof(double));
    
    meanInterpolate(-1, NULL, NUM_POINTS, x, y, TERMS, arr);
	//interpolate(TERMS, x, y, arr);
    printf("\nThe  polynomial is:\n\ny = ");
    for(j = 0; j < TERMS; j++)
    {
        printf(" %+2.4E x^%d ",  arr[j], j);

    }
	printf("\n\n");

    //printf("\n%d\n", binomialCoef(NUM_POINTS, TERMS)	);

    free(x);
    free(y);
    free(arr);

    return 0;
}

