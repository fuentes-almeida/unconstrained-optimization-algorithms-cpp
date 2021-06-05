#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "Configuration.h"
#include "PGMFunctions.h"
#include "ConjGrad.h"
#include "GSeidel.h"

using namespace std;

int main(int argc, char**argv)
{
    char* filename= argv[1];
    double lambda=atof(argv[2]);
    clock_t start,finish;

    Matrix *ImageY=(Matrix*)malloc(sizeof(Matrix));
    readPGM(filename,ImageY);
    Matrix *ImageX1=RequestMatrixDoubleMem(ImageY->rows,ImageY->cols);
    Matrix *ImageX2=RequestMatrixDoubleMem(ImageY->rows,ImageY->cols);

    printf("Prueba con lambda=%lf\n\n",lambda);

    printf("Gradiente Conjugado:\n");
    start=clock();
    ImageConjugatedGradient(ImageX1,ImageY,lambda);
    finish=clock();
    printf("%lf seg\n",(finish-start)/(double)CLOCKS_PER_SEC);
    ///La funcion que calcula el error entre las dos imagenes se encuentra en PGMFunctions
    printf("Error = %lf\n\n",ErrorFunction(ImageX1,ImageY,lambda));

    printf("Gauss Seidel:\n");
    start=clock();
    ImageGaussSeidel(ImageX2,ImageY,lambda);
    finish=clock();
    printf("%lf seg\n",(finish-start)/(double)CLOCKS_PER_SEC);
    ///La funcion que calcula el error entre las dos imagenes se encuentra en PGMFunctions
    printf("Error = %lf\n\n",ErrorFunction(ImageX2,ImageY,lambda));

    FreeMatrixMem(ImageX1);
    FreeMatrixMem(ImageX2);
    FreeMatrixMem(ImageY);

    return 0;
}

