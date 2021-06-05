#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "MatrixInterface.h"
#include "PGMFunctions.h"
#include "dogleg.h"
#include "mrc.h"
#include "functions.h"
using namespace std;


int main(int arc, char**argv)
{
    char *inputfile1=argv[1];
    char *inputfile2=argv[2];

    Matrix *Image1=RequestMatrixDoubleMem(imgsize,imgsize);
    Matrix *ImageRef=RequestMatrixDoubleMem(imgsize,imgsize);
    readPGM(inputfile1,Image1);
    readPGM(inputfile2,ImageRef);

    Matrix **ImgGrad=(Matrix**)malloc(2*sizeof(Matrix*));
    for (int i=0;i<2;i++)
        ImgGrad[i]=RequestMatrixDoubleMem(imgsize,imgsize);
    ImageGradient(Image1,ImgGrad);

    Matrix **ImgHessian=(Matrix**)malloc(4*sizeof(Matrix*));
    for (int i=0;i<4;i++)
        ImgHessian[i]=RequestMatrixDoubleMem(imgsize,imgsize);
    ImageHessian(ImgGrad,ImgHessian);

    TrustRegion(Image1,ImgGrad,ImgHessian,ImageRef);

    FreeMatrixMem(Image1);
    FreeMatrixMem(ImageRef);

    for (int i=0;i<2;i++)
        FreeMatrixMem(ImgGrad[i]);
    free(ImgGrad);

    for (int i=0;i<4;i++)
        FreeMatrixMem(ImgHessian[i]);
    free(ImgHessian);

    return 0;
}
