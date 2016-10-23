// **** Include libraries here ****
// Standard libraries
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//CMPE13 Support Library


// User libraries
#include "../hed/MathMethods.h"

// **** Set macros and preprocessor directives ****

#define  ROWS 4
#define  COLS 4

// **** Define global, module-level, or external variables here ****

// **** Declare function prototypes ****

int main(void)
{
    /******************************************************************************
     * Your code goes in between this comment and the following one with asterisks.
     *****************************************************************************/


    int i;
    int j;
    
    double *x = calloc(5, sizeof(double));
    
    x[0] = 1;
    x[1] = 3;
    x[2] = 4;
    x[3] = -5;
    x[4] = -20;
    //double z[5] = {0, 0, 0, 0, 0};
    
    double *y = calloc(5, sizeof(double));
    
    y[0] = 0.01;
    y[1] = 21.93;
    y[2] = 39.5;
    y[3] = 29.82;
    y[4] = 735.3;

    

    double *arr = calloc(3, sizeof(double));

    //double myMatrix = topRow;
    double myMatrix[ROWS][COLS] = {
        {7, 3, 1, 4},
        {15, 7, 3, 9},
        {6, 5, 3, 1},
        {10, 2, 0, 8},
    };

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

    printf("\n\nI will interpolate these (x, y) pairs into a parabola:\n");
    
    for(j = 0; j < 5; j++)
    {   
        (j<4) ? printf(" ") : 0;  //print an extra space to make it pretty
        printf(" (%+2.4f   ,  %+2.4f)\n", x[j], y[j]);
    }
    
    
    interpolate(5, x, y, 2, arr);
    printf("\nThe  parabola is:\n\n    y = ");
    for(j = 0; j < 3; j++)
    {
        printf(" %2.4f x^(%d) ", (float) arr[j], j);
        if(j < 2)
        {
            printf("+");
        }
    }
    printf("\n\nYou did it, Hooray!\n");

    free(x);
    free(y);
    free(arr);
    return 0;
}
