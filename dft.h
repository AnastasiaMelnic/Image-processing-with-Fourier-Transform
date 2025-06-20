#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define M 8 //inaltime (numar linii)
#define N 8  //latime (numar coloane)
#define PI 3.14159265

typedef unsigned char uint8;
typedef signed char sint8;

//Tablouri pentru imaginea de intrare, imaginile rezultate si spectrele
uint8 imgInput[M*N];  
uint8 imgReconstituitaDirect[M*N];  
uint8 imgOut_Butterworth_LPF[M*N];
uint8 imgOut_Butterworth_HPF[M*N];
uint8 Magnitude_spectrum[M*N];
uint8 Phase_spectrum[M*N];

//Declaratiile functiilor
uint8 round_positive_float(float num);
void DFT_2D(uint8 *imgInput, float *parte_reala_DFT, float* parte_imag_DFT);
void IDFT_2D(float *parte_reala_DFT, float *parte_imag_DFT, uint8 *imgOutput);
void Spectru_modul(float *parte_reala_DFT, float* parte_imag_DFT, float* modul_DFT);
void Spectru_faza(float *parte_reala_DFT, float* parte_imag_DFT, float* faza_DFT);
void Butterworth_LPF(uint8 *imgOutput, float* parte_reala_DFT, float* parte_imag_DFT, float W0);
void Butterworth_HPF(uint8 *imgOutput, float* parte_reala_DFT, float* parte_imag_DFT, float W0);
