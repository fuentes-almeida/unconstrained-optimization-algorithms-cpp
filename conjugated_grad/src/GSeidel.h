#ifndef GSEIDEL_H_INCLUDED
#define GSEIDEL_H_INCLUDED
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "Configuration.h"
#include "PGMFunctions.h"

///Funcion que suaviza la imagen por el metodo de Gauss Seidel
///ImageX=imagen suavizada, ImageY=imagen observada, lambda=parametro de regularizacion
void ImageGaussSeidel(Matrix *ImageX, Matrix *ImageY, double lambda)
{
    double error=1e9,epsilon=1e-6,result,neighbours;
    int num,k=0,Itmax=2*ImageX->size;
    int N=ImageX->rows,M=ImageX->cols;

    const char *outputfile="GSeidel.pgm";
    Matrix *Aux=RequestMatrixDoubleMem(N,M);

    ///valor inicial de la imagen
    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
        {
            num=NumNeighbours(ImageX,i,j);
            ImageX->data[i*M+j]=ImageY->data[i*M+j]/(1+lambda*num);
        }

    while(error>epsilon && k<Itmax)
    {
        k++;error=0;
        for (int i=0;i<N;i++)
            for (int j=0;j<M;j++)
            {
                ///calculo de la formula ya establecida para el metodo Gauss Seidel
                num=NumNeighbours(ImageX,i,j);
                neighbours=PixelNeighbours(ImageX,i,j);

                ///la imagen nueva depende de los vecinos de su iteracion anterior
                result=(ImageY->data[i*M+j]+lambda*neighbours)/(1+lambda*num);
                error+=pow(ImageX->data[i*M+j]-result,2);
                Aux->data[i*M+j]=result;
            }
        Matrix_copy(Aux,ImageX);
        error=sqrt(error);
    }
    printf("%d iteraciones\n",k);

    NormalizeImage(ImageX,ImageX);
    writePGM(outputfile,ImageX);

    FreeMatrixMem(Aux);
}


#endif // GSEIDEL_H_INCLUDED
