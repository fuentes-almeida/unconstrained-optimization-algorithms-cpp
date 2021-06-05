#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include "PGMFunctions.h"
#include "MatrixInterface.h"
#include "CommonTypes.h"
#include "GaussNewton.h"
#include "LevMaquardt.h"
#include "functions.h"

using namespace std;

int main(int argc, char** argv)
{
    srand(time(NULL));
    char *filenamme=argv[1];
    int methodID=atoi(argv[2]);

    Matrix *ImgWrapped=(Matrix*)malloc(sizeof(Matrix));
    readPGM(filenamme,ImgWrapped);
    PiNormalizeImage(ImgWrapped,ImgWrapped);

    NonLinearLeastSquare(ImgWrapped,methodID);

    FreeMatrixMem(ImgWrapped);

    return 0;
}
