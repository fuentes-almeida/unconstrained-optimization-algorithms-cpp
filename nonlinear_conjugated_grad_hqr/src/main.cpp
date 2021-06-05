#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "Configuration.h"
#include "PGMFunctions.h"
#include "GCNL.h"
#include "function.h"
#include "lineSearch.h"

int main(int argc, char**argv)
{
    srand(time(NULL));
    char* filename=argv[1];
    double lambda=atof(argv[2]);
    double k=atof(argv[3]);
    int MD=atoi(argv[4]);

    Matrix *ImageY=(Matrix*)malloc(sizeof(Matrix));
    readPGM(filename,ImageY);
    Matrix *ImageX=RequestMatrixDoubleMem(ImageY->rows,ImageY->cols);

    NonLinearConjugatedGradient(ImageX,ImageY,lambda,k,MD);

    FreeMatrixMem(ImageX);
    FreeMatrixMem(ImageY);

    return 0;
}
