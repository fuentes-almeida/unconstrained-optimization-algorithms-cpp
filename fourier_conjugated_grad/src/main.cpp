#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <complex.h>
#include <fftw3.h>
#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "Configuration.h"
#include "PGMFunctions.h"
#include "ConjGrad.h"
#include "fourier.h"
#include "ConjGradPrec.h"

int main(int argc, char**argv)
{
    char* filename= argv[1];
    double lambda=atof(argv[2]);
    int flag=atoi(argv[3]);

    clock_t start,finish;

    Matrix *ImageY=(Matrix*)malloc(sizeof(Matrix));
    readPGM(filename,ImageY);
    Matrix *ImageX=RequestMatrixDoubleMem(ImageY->rows,ImageY->cols);

    printf("Prueba con lambda=%lf\n",lambda);

    if (flag==1)
    {
        printf("Gradiente Conjugado:\n");
        start=clock();
        ImageConjugatedGradient(ImageX,ImageY,lambda);
        finish=clock();
        printf("%lf seg\n",(finish-start)/(double)CLOCKS_PER_SEC);
        ///La funcion que calcula el error entre las dos imagenes se encuentra en PGMFunctions
        printf("Error = %lf\n",ErrorFunction(ImageX,ImageY,lambda));
    }
    else if (flag==2)
    {
        printf("Precondicionador de Fourier:\n");
        start=clock();
        ImageConjugatedGradientFourier(ImageX,ImageY,lambda);
        finish=clock();
        printf("%lf seg\n",(finish-start)/(double)CLOCKS_PER_SEC);
        ///La funcion que calcula el error entre las dos imagenes se encuentra en PGMFunctions
        printf("Error = %lf\n",ErrorFunction(ImageX,ImageY,lambda));
    }
    else
    {
        printf("Error de Argumento");
    }

    FreeMatrixMem(ImageX);
    FreeMatrixMem(ImageY);

    return 0;
}
