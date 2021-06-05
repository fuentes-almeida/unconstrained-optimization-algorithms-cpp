#ifndef CONJGRADPREC_H_INCLUDED
#define CONJGRADPREC_H_INCLUDED

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
void ImageConjugatedGradientFourier(Matrix *ImageX, Matrix *ImageY, double lambda)
{
    int k=0, N=ImageX->rows, M=ImageX->cols;
    const char *outputfile="FourierConjGrad.pgm";

    Matrix *gk=RequestMatrixDoubleMem(N,M);
    Matrix *yk=RequestMatrixDoubleMem(N,M);
    Matrix *yk_New=RequestMatrixDoubleMem(N,M);
    Matrix *gk_New=RequestMatrixDoubleMem(N,M);
    Matrix *pk=RequestMatrixDoubleMem(N,M);
    Matrix *aux=RequestMatrixDoubleMem(N,M);

    ///inicializacion de x
    matrixFill(ImageX,50.0);
    ///producto Ax0
    ImageVectorProduct(ImageX,gk,lambda);
    ///se inicia con g0=Ax0-b
    matrixSubstract(gk,ImageY,gk);
    ///Se resuelve para y0
    FourierPrecond(yk->data,gk->data,N,M,lambda);
    ///se calcula p0=-y0
    matrixScalarProduct(pk,yk,-1);

    double epsilon=1e-6,alpha=0,beta=0;
    int Itmax=2*ImageX->size;

    while ((matrixNorm(gk)>epsilon) && (k<Itmax))
    {
        k++;
        ///calculo de alpha
        double gkTyk=matrixInnerProduct(gk,yk);
        double pkTApk=ImageQuadraticFormProduct(pk,lambda);
        alpha=gkTyk/pkTApk;

        ///calculo de la nueva iteracion de x
        matrixScalarProduct(aux,pk,alpha);
        matrixAdd(ImageX,aux,ImageX);

        ///calculo del nuevo gradiente
        matrixFill(aux,0);
        ImageVectorProduct(pk,aux,lambda);
        matrixScalarProduct(aux,aux,alpha);
        matrixAdd(gk,aux,gk_New);

        ///calculo de la nueva yk
        FourierPrecond(yk_New->data,gk_New->data,N,M,lambda);

        ///calculo de beta
        double gkTyk_New=matrixInnerProduct(gk_New,yk_New);
        beta=gkTyk_New/gkTyk;

        ///calculo de la nueva pk
        matrixScalarProduct(pk,pk,beta);
        matrixSubstract(pk,yk_New,pk);

        Matrix_copy(gk_New,gk);
        Matrix_copy(yk_New,yk);
    }

    FreeMatrixMem(gk);
    FreeMatrixMem(gk_New);
    FreeMatrixMem(yk);
    FreeMatrixMem(yk_New);
    FreeMatrixMem(pk);
    FreeMatrixMem(aux);

    printf("%d iteraciones\n",k);

    writePGM(outputfile,ImageX);
}

#endif // CONJGRADPREC_H_INCLUDED
