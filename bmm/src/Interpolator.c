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

#define INTERPOLATOR_FILE_LOCATION "fio/InterpolatorPoints.txt"
#define FILE_MAX_LINE_LEN 1024
#define MAX_NUM_POINTS 512
#define X_ARRAY 1
#define Y_ARRAY 2
#define NO_DATA 3
// **** Declare function prototypes ****
int fscanLine(FILE *filename, char *str)
{
    char currentChar = '\0';
    int i = 0;
    
    while(fscanf(filename, "%c", &currentChar) != EOF)
    {
        
        //if(fscanf(filename, "%c", &currentChar) == EOF)
        //{
          //  return EOF;
        //}
        str[i] = currentChar;
        i++;
        if(currentChar == '\n')
        {
            str[i] = '\0';
            return SUCCESS;
        }
    }
    return EOF;  
   
}


int main(int argc, char **argv)
{
    
    //flush(stdout);
    FILE* pointsFile = (FILE*) fopen( INTERPOLATOR_FILE_LOCATION, "r" );
    
    
    if(pointsFile == NULL)
    {
        printf("\nPlease create file \"%s\" to run this program.\n", 
                INTERPOLATOR_FILE_LOCATION);
        return 1;
    }
    
    int i = 0;
    int numPoints = 0;
    int terms = 0;
    double *x = (double*) calloc(MAX_NUM_POINTS, sizeof(double));
    double *y = (double*) calloc(MAX_NUM_POINTS, sizeof(double));    
    double *arr = (double*) calloc(MAX_NUM_POINTS, sizeof(double));
    char *fileLine = (char*) calloc(FILE_MAX_LINE_LEN , sizeof(char));
    char *tokenizer;
    
    while(fscanLine(pointsFile, fileLine) != EOF)
    { 
            switch(fileLine[0])
            {
                case '~' :
                case '\n' :
                    //tokenizeFlag = NO_DATA;
                    break;
                case 't' :
                    tokenizer = strtok(fileLine, "t: \n");
                    terms = strtod(tokenizer, NULL);
                    break;
                case 'x' :
                    //tokenizeFlag = X_ARRAY;
                    if(fscanLine(pointsFile, fileLine) != EOF)
                    {
                        if(fileLine[strlen(fileLine) - 2] == ' ')
                        {
                            printf("\nPlease format file \"%s\" properly to run" 
                               " this program. Don't leave spaces at the end of"
                               " the list of x values.\n", 
                              INTERPOLATOR_FILE_LOCATION);
                            return 1;
                        }
                        tokenizer = strtok(fileLine, " ");
                        i = 0;
                        while(tokenizer != NULL)
                        {   
                            x[i] = strtod(tokenizer, NULL);
                            tokenizer = strtok(NULL, " ");
                            i++;
                        }
                        numPoints = i;
                    }
                    break;
                case 'y' :
                    //tokenizeFlag = Y_ARRAY;
                    if(numPoints == 0)
                    {
                        printf("\nPlease format file \"%s\" properly to run" 
                               " this program. Refer to README.md for more"
                               " details.\n", 
                              INTERPOLATOR_FILE_LOCATION);
                        return 1;
                    }
                    if(fscanLine(pointsFile, fileLine) != EOF)
                    {
                        if(fileLine[strlen(fileLine) - 2] == ' ')
                        {
                            printf("\nPlease format file \"%s\" properly to run" 
                               " this program. Don't leave spaces at the end of"
                               " the list of y values.\n", 
                              INTERPOLATOR_FILE_LOCATION);
                            return 1;
                        }
                        tokenizer = strtok(fileLine, " ");
                        i = 0;
                        while(tokenizer != NULL && i < numPoints)
                        {
                            y[i] = strtod(tokenizer, NULL);
                            tokenizer = strtok(NULL, " ");
                            i++;
                        }
                        if(i < numPoints)
                        {
                            numPoints = i;
                        }
                    }
                    break;
                default :
                    printf("\nPlease format file \"%s\" properly to run this"
                            " program. Refer to README.md for more details.\n", 
                              INTERPOLATOR_FILE_LOCATION);
                    return 1;
            }
    }
    printf("\nTerms: %d", terms);
    printf("\nX values: ");
    for(i = 0; i < numPoints; i++)
    {
        printf("%2.0f ", x[i]);
    }
    printf("\nY values: ");
    for(i = 0; i < numPoints; i++)
    {
        printf("%2.0f ", y[i]);
    }
    printf("\n\n");

    meanPolynomial(-1, NULL, numPoints, x, y, terms, arr);
	//interpolate(terms, x, y, arr);
    printf("\nThe  polynomial is:\n\ny = ");
    for(i = 0; i < terms; i++)
    {
        printf(" %+2.4E x^%d ",  arr[i], i);

    }
	printf("\n\n");
    

    fclose(pointsFile);
    free(fileLine);
    free(x);
    free(y);    
    free(arr);
    return 0;
}
