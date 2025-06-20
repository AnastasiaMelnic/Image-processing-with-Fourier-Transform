#include "dft.h"

void main()
{


    float parte_reala_DFT[M*N];
    float parte_imag_DFT[M*N];
    float modul_DFT[M*N];
    float faza_DFT[M*N];

	DFT_2D(imgInput, parte_reala_DFT, parte_imag_DFT);

    Spectru_modul(parte_reala_DFT, parte_imag_DFT, modul_DFT);
	Spectru_faza(parte_reala_DFT, parte_imag_DFT, faza_DFT);

	IDFT_2D(parte_reala_DFT, parte_imag_DFT, imgReconstituitaDirect);
	
	Butterworth_LPF(imgOut_Butterworth_LPF, parte_reala_DFT, parte_imag_DFT, 1);
	Butterworth_HPF(imgOut_Butterworth_HPF, parte_reala_DFT,parte_imag_DFT, 1);


}

//Functie ce rotunjeste un numar pozitiv de tip float
uint8 round_positive_float(float num) {
    float int_part = (uint8)num;
    float frac_part = num - int_part;

    if (frac_part >= 0.5f){
        return (uint8)(int_part + 1);
    }
    else{
        return (uint8)int_part;
    }
}

void DFT_2D(uint8 *imgInput, float *parte_reala_DFT, float* parte_imag_DFT) {
	float real_suma, imag_suma, modul, faza;
	int u,v,x,y;
	float imgInput_translatata[M*N];

	sint8 factor_translatare;

	//Translatarea imaginii initiale, pentru a centra spectrul si a putea aplica FTJ/FTS
	for (x = 0; x < M; x++) {
         for (y = 0; y < N; y++) {
             factor_translatare=((x+y)%2==0)?1:-1;
             imgInput_translatata[x * N + y]=(float)imgInput[x * N + y]*factor_translatare;
         }
	}

    //Calcul DFT
    for (u = 0; u < M; u++) {
        for (v = 0; v < N; v++) {
            real_suma = 0.0f;
            imag_suma = 0.0f;

            for (x = 0; x < M; x++) {
                for (y = 0; y < N; y++) {
                    modul = imgInput_translatata[x * N + y];
                    faza = -2.0f * PI * ((float)(u * x) / M + (float)(v * y) / N);
                    real_suma += modul * cosf(faza);
                    imag_suma += modul * sinf(faza);
                }
            }

            parte_reala_DFT[u * N + v] = real_suma;
            parte_imag_DFT[u * N + v] = imag_suma;
        }
    }


}

//Functie ce calculeaza inversa Transformatei Fourier discrete, care realizeaza trecerea din 
//domeniul spectral in domeniul spatial
void IDFT_2D(float *parte_reala_DFT, float *parte_imag_DFT, uint8 *imgOutput) {
    int x, y, u, v;
    float real_suma;
    sint8 factor_translatare;

    float img_temp[M * N];
	//Calcul IDFT
    for (x = 0; x < M; x++) {
        for (y = 0; y < N; y++) {
            real_suma = 0.0f;

            for (u = 0; u < M; u++) {
                for (v = 0; v < N; v++) {
                    float faza = 2.0f * PI * ((float)(u * x) / M + (float)(v * y) / N);
                    float real = parte_reala_DFT[u * N + v];
                    float imag = parte_imag_DFT[u * N + v];

                    real_suma += real * cosf(faza) - imag * sinf(faza);
                }
            }
			//Translatare
            factor_translatare = ((x + y) % 2 == 0) ? 1 : -1;
            img_temp[x * N + y] = real_suma * factor_translatare/(M*N);
        }
    }

	
    for (x = 0; x < M; x++) {
        for (y = 0; y < N; y++) {
            float valoare = img_temp[x * N + y];
            //Pentru a fi siguri ca nu exista valori in afara intervalui [0 255]
            if (valoare < 0.0f) valoare = 0.0f;
            if (valoare > 255.0f) valoare = 255.0f;
            //Copiaza valorile din tabloul img_temp[] in imgOutput[]
            imgOutput[x * N + y] = (uint8)round_positive_float(valoare);
        }
	}
}

//Calculeaza spectrul de amplitudine
void Spectru_modul(float *parte_reala_DFT, float* parte_imag_DFT, float* modul_DFT){
	int i;
	float max_modul=0.0f;
	for(i=0; i<M*N; i++){
		modul_DFT[i]=sqrtf(parte_reala_DFT[i]*parte_reala_DFT[i]+parte_imag_DFT[i]*parte_imag_DFT[i]);
		modul_DFT[i]=logf(1+modul_DFT[i]);
		if(modul_DFT[i]>max_modul){
			max_modul=modul_DFT[i];
		}
	}
	//Normalizare valori pentru afisare
	for(i=0; i<M*N; i++){
		modul_DFT[i]/=max_modul;
		modul_DFT[i]*=255.0f;
		//Functia round_positive_float rotunjeste valorile, pentru a nu fi trunchiate
		Magnitude_spectrum[i]=(uint8)round_positive_float(modul_DFT[i]);
	}
}
//Calculeaza spectrul de faza
void Spectru_faza(float *parte_reala_DFT, float* parte_imag_DFT, float* faza_DFT){
	int i;
	float faza_min,faza_max;

	faza_DFT[0]=atan2f(parte_imag_DFT[0], parte_reala_DFT[0]);
	faza_min = faza_DFT[0];
	faza_max = faza_DFT[0];

	//Normalizare valori pentru afisare
	for(i=1; i<M*N; i++){
		faza_DFT[i]=atan2f(parte_imag_DFT[i], parte_reala_DFT[i]);
		if (faza_DFT[i] < faza_min){
			 faza_min = faza_DFT[i];
		}
        if (faza_DFT[i] > faza_max){
        	 faza_max = faza_DFT[i];
        }
	}
	for (i = 0; i < M * N; i++) {
        faza_DFT[i] = 255.0f * (faza_DFT[i] - faza_min) / (faza_max - faza_min);
		Phase_spectrum[i]=(uint8)round_positive_float(faza_DFT[i]);
	}
}

//Aplica filtru trece jos Butterworth
void Butterworth_LPF(uint8 *imgOutput, float* parte_reala_DFT, float* parte_imag_DFT, float W0) {
    int u, v;
    float W, H, du, dv, parte_reala_filtrata[M*N], parte_imag_filtrata[M*N];


    for (u = 0; u < M; u++) {
        for (v = 0; v < N; v++) {

            float du = (float)(u - M / 2.0f);
            float dv = (float)(v - N / 2.0f);
            W = sqrtf((float)(du * du + dv * dv));


            H = 1.0f / (1.0f + powf((W / W0),2.0f));

            parte_reala_filtrata[u * N + v]=parte_reala_DFT[u * N + v] * H;
            parte_imag_filtrata[u * N + v]=parte_imag_DFT[u * N + v] * H;
        }
    }
    //Calculeaza inversa Transformatei Fourier Discrete pentru a putea vizualiza imaginea rezultata 
    //in urma filtrarii
	IDFT_2D(parte_reala_filtrata, parte_imag_filtrata,imgOutput);
}

//Aplica filtru trece sus Butterworth
void Butterworth_HPF(uint8 *imgOutput, float* parte_reala_DFT, float* parte_imag_DFT, float W0) {
    int u, v;
    float W, H, du, dv, parte_reala_filtrata[M*N], parte_imag_filtrata[M*N];


    for (u = 0; u < M; u++) {
        for (v = 0; v < N; v++) {

            float du = (float)(u - M / 2.0f);
            float dv = (float)(v - N / 2.0f);
            W = sqrtf((float)(du * du + dv * dv));

            if(W==0){
            	H=0.0f;
            }
            else{
                H = 1.0f / (1.0f + powf((W0 / W),2.0f));
            }

            parte_reala_filtrata[u * N + v]=parte_reala_DFT[u * N + v] * H;
            parte_imag_filtrata[u * N + v]=parte_imag_DFT[u * N + v] * H;
        }
    }
    //Calculeaza inversa Transformatei Fourier Discrete pentru a putea vizualiza imaginea rezultata 
    //in urma filtrarii
	IDFT_2D(parte_reala_filtrata, parte_imag_filtrata,imgOutput);
}

