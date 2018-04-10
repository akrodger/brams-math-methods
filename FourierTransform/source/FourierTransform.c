
//GBA Libs

#include <gba_console.h>
#include <gba_video.h>
#include <gba_interrupt.h>
#include <gba_systemcalls.h>
#include <gba_input.h>
#include <mbv2.h>
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
#define NUM_FREQS  5



void VBLANK_ISR_FUNCTION() {
    scanKeys();
    //REG_IF = IRQ_VBLANK;
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

	//sineWave: the function to be approximated
	double *sineWave = calloc(SCREEN_WIDTH, sizeof(double));
	//xCoords: the x coordinates of the graphs.
	double *xCoords = calloc(SCREEN_WIDTH, sizeof(double));
	//fourierCos: the coeficients of cosines in the Fourier Transform
	double *fourierCos = calloc(NUM_FREQS, sizeof(double));
	//fourierSin: the coeficients of sines in the Fourier Transform
	double *fourierSin = calloc(NUM_FREQS, sizeof(double));
	//amps: amplitudes of each term (randomized each loop)
	double *amps = calloc(NUM_FREQS*2, sizeof(double));
	//iterators
    int i = 0;
	int j = 0;
	unsigned int randCounter = 0;
	double waveValue = 0;
	double errorValue = 0;
	double meanOfData = 0;
	consoleDemoInit();
	iprintf("\x1b[6;1HThis demo generates");
	iprintf("\x1b[7;1Hnoisy waveforms (red),");
	iprintf("\x1b[8;1Hthen computes their");
	iprintf("\x1b[9;1HFourier Projections (blue).");
	iprintf("\x1b[12;4HPress Any Key.");
	while(keysDown()==0){
		randCounter++;
	}
	srand(randCounter);
	iprintf("\x1b[14;8HWorking...");
	while(1){

		//generate amplitudes -1 to 1
		for(i = 0; i < NUM_FREQS*2; i++){
			amps[i] = ((double)rand()/RAND_MAX*2.0-1.0);
		}
		meanOfData = ((double)rand()/RAND_MAX*2.0-1.0);
		//evaluate all the trig functions at the x coordinates
		//store those value in sineWave
		for(i = 0; i < SCREEN_WIDTH; i++){
			sineWave[i] = 0;
			for(j = 0; j < NUM_FREQS; j++){
				sineWave[i] += 
						amps[2*j]*(cos(M_PI * 2*(j+1)*i / SCREEN_WIDTH));
				sineWave[i] +=
						amps[(2*j)+1]*(sin(M_PI*2*(j+1)*i/SCREEN_WIDTH));
			}
			//make some NOISE! WOOOOO
			errorValue = ((double)rand()/RAND_MAX*0.5);
			sineWave[i] += errorValue*((double)rand()/RAND_MAX*2.0-1.0);

			sineWave[i] += meanOfData;

			xCoords[i] = ((double) i)/(SCREEN_WIDTH);
		}
		meanOfData = discreteFourier(SCREEN_WIDTH,
						xCoords,
						sineWave,
						NUM_FREQS,
						fourierCos,
						fourierSin);
		//set display to background drawing mode.
		if(REG_DISPCNT != (BG2_ON | MODE_3)){
			REG_DISPCNT = BG2_ON | MODE_3;
		}

		for(i = 0; i < SCREEN_WIDTH; i++){
			for(j = 0; j < SCREEN_HEIGHT; j++){
			LCD_WRITE(i, j, 0, 0, 0);
			}
		}
		for(i = 0; i < SCREEN_WIDTH; i++){
			LCD_WRITE(i, (int) (((-SCREEN_HEIGHT/8)*sineWave[i]) +
								(SCREEN_HEIGHT/2)), 255, 0, 0);
			waveValue = 0;
			for(j = 0; j < NUM_FREQS; j++){
				waveValue +=
				(sin(M_PI * 2*(j+1)*i / SCREEN_WIDTH))*fourierSin[j]+
				(cos(M_PI * 2*(j+1) * i / SCREEN_WIDTH))*(fourierCos[j]);
			}
			waveValue += meanOfData;
			LCD_WRITE(i,
				(int) (-(SCREEN_HEIGHT/8)*(waveValue) + SCREEN_HEIGHT/2),
						0, 0, 255);
		}
	}
}
