#ifndef CONJGRAD_H_INCLUDED
#define CONJGRAD_H_INCLUDED

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "PGMFunctions.h"

///funcion que calcula el producto matriz-vector, Image es la imagen que se esta suevizando,
///Vec es el resultado del producto y scalar es el valor de lambda
void ImageVectorProduct(Matrix *Image, Matrix *Vec, double scalar)
{
    int M=Image->cols;
    int N=Image->rows;

    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
        {
            ///Las funciones NumNeighbours() y PixelNeighbours() son comunes a los dos metodos,
            ///se encuentran en el archivo de PGMFunctions, calculan el numero de vecinos y la suma, respectivamente
            int num=NumNeighbours(Image,i,j);
            int sum=PixelNeighbours(Image,i,j);
            Vec->data[i*M+j]=(1+scalar*num)*Image->data[i*M+j]-scalar*sum;
        }
}

///calcula el producto cuadratico vector-matriz-vector de forma eficiente
double ImageQuadraticFormProduct(Matrix *Vec, double lambda)
{
    Matrix *Aux=RequestMatrixDoubleMem(Vec->rows,Vec->cols);
    ImageVectorProduct(Vec,Aux,lambda);
    double product=matrixInnerProduct(Aux,Vec);
    FreeMatrixMem(Aux);
	return product;
}

///Funcion que suaviza la imagen por el metodo de Gradiente Conjugado
///ImageX=imagen suavizada, ImageY=imagen observada, lambda=parametro de regularizacion
void ImageConjugatedGradient(Matrix *ImageX, Matrix *ImageY, double lambda)
{
    int k=0;
    const char *outputfile="ConjGrad.pgm";
    Matrix *gradk=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);
    Matrix *gradk_New=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);
    Matrix *pk=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);
    Matrix *aux=RequestMatrixDoubleMem(ImageX->rows,ImageX->cols);

    ImageVectorProduct(ImageX,gradk,lambda);///equivale a matriz x vector
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

    NormalizeImage(ImageX,ImageX);
    writePGM(outputfile,ImageX);
}



#endif // CONJGRAD_H_INCLUDED
