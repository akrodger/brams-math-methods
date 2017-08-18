
//GBA Libs

#include <gba_console.h>
#include <gba_video.h>
#include <gba_interrupt.h>
#include <gba_systemcalls.h>
#include <gba_input.h>	
//Standard Libs
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//User Libs
#include "MathMethods.h"
#include "Stack.h"


typedef void (*fnptr)(void);
#define REG_ISR_MAIN *(fnptr*)(0x03007FFC)

#define LCD_WRITE(x, y, r, g, b) \
                    ((unsigned short*)0x06000000)[(x)+(y)*SCREEN_WIDTH]  = \
                                                            RGB8((r), (g), (b))
#define NUM_FUNS 15

//These are the functions to be used with LeastSquares Algorithm
void VBLANK_ISR_FUNCTION() {
    scanKeys();
    //REG_IF = IRQ_VBLANK;
}

double constFun(double x){
	return 1;
}

double lineFun(double x){
	return x;
}

double quadFun(double x){
	return x*x;
}

double cubeFun(double x){
	return x*x*x;
}

double quarFun(double x){
	return x*x*x*x;
}

double quinFun(double x){
	return x*x*x*x*x;
}

double sectFun(double x){
	return x*x*x*x*x*x;
}

double septFun(double x){
	return x*x*x*x*x*x*x;
}

double octoFun(double x){
	return x*x*x*x*x*x*x*x;
}

double nonaFun(double x){
	return x*x*x*x*x*x*x*x*x;
}

double decaFun(double x){
	return x*x*x*x*x*x*x*x*x*x;
}

double elevFun(double x){
	return x*x*x*x*x*x*x*x*x*x*x;
}
double twelFun(double x){
	return x*x*x*x*x*x*x*x*x*x*x*x;
}
double thirtFun(double x){
	return x*x*x*x*x*x*x*x*x*x*x*x*x;
}

double fourtFun(double x){
	return x*x*x*x*x*x*x*x*x*x*x*x*x;
}
//---------------------------------------------------------------------------------
// Program entry point
//---------------------------------------------------------------------------------
int main(void) {
//---------------------------------------------------------------------------------


	// the vblank interrupt must be enabled for VBlankIntrWait() to work
	// since the default dispatcher handles the bios flags no vblank handler
	// is required
	irqInit();
	irqEnable(IRQ_VBLANK);
    irqSet(IRQ_VBLANK, VBLANK_ISR_FUNCTION);

	//polyWave: The function to be approximated
	double *polyWave = calloc(SCREEN_WIDTH, sizeof(double));
	//xCoords: the x coordinates of the graphs.
	double *xCoords = calloc(SCREEN_WIDTH, sizeof(double));
	//polyCoefs: the approximation coeficients to be solved for
	double *polyCoefs = calloc(NUM_FUNS, sizeof(double));
	//amps: amplitudes of each term (randomized each loop)
	double *amps = calloc(NUM_FUNS, sizeof(double));
	//phi: array of double->double function pointers
	double (*phi[NUM_FUNS])(double);
	//initialize all the elements of phi
	phi[0] = constFun;
	phi[1] = lineFun;
	phi[2] = quadFun;
	phi[3] = cubeFun;
	phi[4] = quarFun;
	phi[5] = quinFun;
	phi[6] = sectFun;
	phi[7] = septFun;
	phi[8] = octoFun;
	phi[9] = nonaFun;
	phi[10] = decaFun;
	phi[11] = elevFun;
	phi[12] = twelFun;
	phi[13] = thirtFun;
	phi[14] = fourtFun;
	//iterators
    int i = 0;
	int j = 0;
	unsigned int randCounter = 0;
	//polyValue: used to sum over the coeficients of a polynomial for graphing 
	double polyValue = 0;
	double errorValue = 0;
	double meanOfData = 0;
	consoleDemoInit();	
	iprintf("\x1b[6;1HThis demo generates a");
	iprintf("\x1b[7;1Hnoisy waveform (red),");
	iprintf("\x1b[8;1Hthen computes a 15 term");
	iprintf("\x1b[9;1Hbest-fit polynomial (blue).");
	iprintf("\x1b[12;4HPress Any Key.");
	while(keysDown()==0){
		randCounter++;
	}	
	srand(randCounter);
	iprintf("\x1b[14;8HWorking...");
	while(1){
		//generate random amplitudes -1 to 1
		for(i = 0; i < NUM_FUNS; i++){
			amps[i] = ((double)rand()/RAND_MAX*1.2-0.6	);
		}
		meanOfData = ((double)rand()/RAND_MAX*2.0-1.0);
		//evaluate a trig series. We approximate the trig series by a polynomial
		//store those values into polyWave
		for(i = 0; i < SCREEN_WIDTH; i++){
			xCoords[i] = (2.0*((double) i)/(SCREEN_WIDTH)) - 1.0;
			polyWave[i] = 0;
			for(j = 0; j < NUM_FUNS/5; j++){
				polyWave[i] += amps[2*j]*sin(2*M_PI*(j)*xCoords[i]);
				polyWave[i] += amps[2*j +1]*cos(2*M_PI*(j)*xCoords[i]);
			
			}  
			//this line adds randomness to polyWave (make some NOISEEEEE)
			errorValue = ((double)rand()/RAND_MAX*0.25);
			polyWave[i] += errorValue*((double)rand()/RAND_MAX*2.0-1.0);
			polyWave[i] += meanOfData;
		}
		//compute an approximation of polyWave, store the coeficients of
		//that approximation into polyCoefs
		//this function takes a few seconds to compute, giving you the ability
		//to view whatever was computed on the previous iteration of the loop.
		discreteLeastSquares(SCREEN_WIDTH,
							xCoords,
							polyWave,
							NUM_FUNS,
							phi, 
							polyCoefs);
		//set display to background drawing mode.
		if(REG_DISPCNT != (BG2_ON | MODE_3)){
			REG_DISPCNT = BG2_ON | MODE_3;
		}
		//Clear the screen
		for(i = 0; i < SCREEN_WIDTH; i++){
			for(j = 0; j < SCREEN_HEIGHT; j++){
				LCD_WRITE(i, j, 0, 0, 0);					
			}
		}
		//Draw the randomized polyWave on screen in red, then draw the 
		//computed approximation with polyCoefs on screen in blue.
		for(i = 0; i < SCREEN_WIDTH; i++){
			LCD_WRITE(i, (int) (((-SCREEN_HEIGHT/4)*polyWave[i]) + 
								(SCREEN_HEIGHT/2)), 255, 0, 0);
			polyValue = 0;
			for(j = 0; j < NUM_FUNS; j++){
				polyValue += polyCoefs[j]*(*phi[j])(xCoords[i]);
			}
			
			LCD_WRITE(i,
					(int) (-(SCREEN_HEIGHT/4)*(polyValue) + SCREEN_HEIGHT/2), 
					0, 0, 255);
		}
		//end of loop. 
	}
}


