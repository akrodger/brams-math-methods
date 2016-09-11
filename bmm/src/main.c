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
            "\nOur Original Matrix. Here we will call it A:\n"
            "This matrix is something I wrong ksdjflsakdg");
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
    /******************************************************************************
     * Your code goes in between this comment and the preceding one with asterisks
     *****************************************************************************/

    return 0;
}
