#ifndef FOURIER_H_INCLUDED
#define FOURIER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex.h>
#include <fftw3.h>
#include "MatrixInterface.h"
#include "CommonTypes.h"

void FourierPrecond(double *Fxy, double *Gxy, int N, int M, double lambda)
{
    fftw_complex *DFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*M);

    fftw_plan direct = fftw_plan_dft_2d(N, M, DFT, DFT, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan inverse = fftw_plan_dft_2d(N, M, DFT, DFT, FFTW_BACKWARD, FFTW_ESTIMATE);

    ///vaciar Gxy en el complejo
    for (int i=0;i<N*M;i++){
        DFT[i][0]=Gxy[i];
        DFT[i][1]=0.0;
    }
    ///Se ejecuta la transformada de Fourier
    fftw_execute(direct);

    ///Se aplica la formula del precondicionador
    for (int u=0; u<N; u++)
        for (int v=0; v<M; v++)
        {
            double Huv=1+4*lambda-2*lambda*(cos(2*PI*u/N)+cos(2*PI*v/M));
            DFT[u*M+v][0]=DFT[u*M+v][0]/Huv;
            DFT[u*M+v][1]=DFT[u*M+v][1]/Huv;
        }

    ///Se ejecuta la transformada inversa de Fourier
    fftw_execute(inverse);

    ///vaciar el resultado en Fxy
    for (int i=0;i<N*M;i++)
        Fxy[i]=DFT[i][0]/(N*M);

    fftw_destroy_plan(direct);
    fftw_destroy_plan(inverse);
    fftw_free(DFT);
}

#endif // FOURIER_H_INCLUDED
