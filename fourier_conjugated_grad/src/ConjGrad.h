#ifndef CONJGRAD_H_INCLUDED
#define CONJGRAD_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fftw3.h>

#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "PGMFunctions.h"

///Funcion que suaviza la imagen por el metodo de Gradiente Conjugado
///ImageX=imagen suavizada, ImageY=imagen observada, lambda=parametro de regularizacion
///Las funciones PixelNeighbours(), ImageVectorProduct(), e ImageQuadraticFormProduct() son comunes a los metodos,
///se encuentran en el archivo de PGMFunctions, calculan la suma de los vecinos y los productos matriz-vector eficientes
void ImageConjugatedGradient(Matrix *ImageX, Matrix *ImageY, double lambda)
{
    int k=0;
    const char *outputfile="ConjGrad.pgm";
    Matrix *gradk=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);
    Matrix *gradk_New=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);
    Matrix *pk=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);
    Matrix *aux=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);

    matrixSubstract(gradk,ImageY,gradk);
    matrixScalarProduct(pk,gradk,-1);

    double epsilon=1e-6,alpha=0,beta=0;
    int Itmax=2*ImageX->size;

    while ((matrixNorm(gradk)>epsilon) && (k<Itmax))
    {
        k++;
        ///calculo de alpha
        double gkTpk=matrixInnerProduct(gradk,pk);
        double pkTApk=ImageQuadraticFormProduct(pk,lambda);
        alpha=-gkTpk/pkTApk;

        ///calculo de la nueva iteracion de x
        matrixScalarProduct(aux,pk,alpha);
        matrixAdd(ImageX,aux,ImageX);

        ///calculo del nuevo gradiente
        matrixFill(aux,0);
        ImageVectorProduct(pk,aux,lambda);
        matrixScalarProduct(aux,aux,alpha);
        matrixAdd(gradk,aux,gradk_New);

        ///calculo de beta
        double gkTgk_New=matrixInnerProduct(gradk_New,gradk_New);
        double gkTgk=matrixInnerProduct(gradk,gradk);
        beta=gkTgk_New/gkTgk;

        ///calculo de la nueva pk
        matrixScalarProduct(pk,pk,beta);
        matrixSubstract(pk,gradk_New,pk);

        Matrix_copy(gradk_New,gradk);
    }

    FreeMatrixMem(gradk);
    FreeMatrixMem(gradk_New);
    FreeMatrixMem(pk);
    FreeMatrixMem(aux);

    printf("%d iteraciones\n",k);

    writePGM(outputfile,ImageX);
}

#endif // CONJGRAD_H_INCLUDED
