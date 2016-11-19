// **** Include libraries here ****
// Standard libraries
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


// User libraries
#include "MathMethods.h"
#include "Stack.h"



int main(int argc, char **argv)
{

    int i;
    int j;

    if(argc < 3 || 4 < argc)
	{
		printf("\nUsage: LinAlg [ROWS] [COLS] [BUFFER]\n"
		"NOTE: 3rd argument is optional floating point percision setting.\n\n");
		return 1;
	}

    char *end = NULL;
	int ROWS = strtol(argv[1], &end, 0);
    
    int COLS = strtol(argv[2], &end, 0);

    if(argv[3] != NULL)
	{
		double deltaBuffer = strtod(argv[3], &end);
		setDelta(deltaBuffer);
	}

    //double myMatrix = topRow;
    double myMatrix[ROWS][COLS];
    
    for(i = 0; i < ROWS; i++)
    {
        for(j = 0; j < COLS; j++)
        {
            myMatrix[i][j] = pow(i, j);
        }
    }

    double copiedMat[ROWS][COLS];

    double product[ROWS][COLS];


    copy(ROWS, COLS, myMatrix, copiedMat);


    double det = ref(ROWS, COLS, copiedMat);


    printf("\nThis main serves a short testing file.\n"
            "An example of the library is shown here.\n"
            "\nOur Original Matrix. Here we will call it A:\n");
    for (i = 0; i < ROWS; i++) {
        printf("\n");
        for (j = 0; j < COLS; j++) {
            //*(topRow + ((int) sizeof(double)) * i + j) = ;
            printf(" %+2.4f ", (float) myMatrix[i][j]);
        }
    }

    printf("\n\nUsing ref(), this is a Row Echelon Form of A.\n");
    for (i = 0; i < ROWS; i++) {
        printf("\n");
        for (j = 0; j < COLS; j++) {
            printf(" %+2.4f ", (float) copiedMat[i][j]);
        }
    }

    printf("\n\nAnd this is its determinant %+2.4f according to ref().",
            (float) det);

    copy(ROWS, COLS, myMatrix, copiedMat);

    det = invert(ROWS, COLS, copiedMat);

    printf("\n\nThis is the inverse of A:\n");
    for (i = 0; i < ROWS; i++) {
        printf("\n");
        for (j = 0; j < COLS; j++) {
            printf(" %+2.4f ", (float) copiedMat[i][j]);
        }
    }

    printf("\n\nIts determinant is as calculated by invert() is %+2.4f.",
            (float) det);

    printf("\n\nThis is the Product of A and its inverse:\n");

    multiply(ROWS, COLS, COLS, myMatrix, copiedMat, product);

    for (i = 0; i < ROWS; i++) {
        printf("\n");
        for (j = 0; j < COLS; j++) {
            printf(" %+2.4f ", (float) product[i][j]);
        }
    }


    printf("\n\nYou did it, Hooray!\n");

    return 0;
}
