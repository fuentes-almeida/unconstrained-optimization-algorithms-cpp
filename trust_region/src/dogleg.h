#ifndef DOGLEG_H_INCLUDED
#define DOGLEG_H_INCLUDED

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "MatrixInterface.h"
#include "PGMFunctions.h"
#include "functions.h"
using namespace std;

///calcula pk por el metodo de dogleg
void Dogleg(double delta, Vector *gradient, Matrix *Hess, Vector *pC)
{
    vectorFill(pC,0);
    Vector *gTB=RequestVectorDoubleMem(dim);
    Vector *pU=RequestVectorDoubleMem(dim);
    Vector *pB=RequestVectorDoubleMem(dim);
    Vector *aux=RequestVectorDoubleMem(dim);
    Matrix *Bk=RequestMatrixDoubleMem(dim,dim);
    Matrix_copy(Hess,Bk);
    double lambda=10;

    while (quadraticFormProduct(Bk,gradient)<0)
        matrixAddToDiagonal(Bk,Bk,lambda);

    ///calculo de pU: de Gradiente
    double gTg=vectorInnerProduct(gradient,gradient);
    double gTBg=quadraticFormProduct(Bk,gradient);
    vectorScalarProduct(pU,gradient,-(gTg/gTBg));
    double pUNorm=vectorNorm(pU);

    ///calculo de pB: de Newton
    vectorScalarProduct(aux,gradient,-1);
    SolveCrouthSystem(Bk,aux,pB,dim);

    ///calculo de pC: punto de Cauchy:

    if (vectorNorm(pB)<=delta)
        {
            vectorCopy(pB,pC);
            cout<<"pB: ";
        }

    else if (vectorNorm(pU)>=delta)
    {
        vectorScalarProduct(pC,pU,delta/pUNorm);
        cout<<"pU ";
    }


    else
    {
        Vector *pB_pU=RequestVectorDoubleMem(dim);
        vectorSubstract(pB,pU,pB_pU);

        double a=vectorNorm2(pB_pU);
        double b=2*vectorInnerProduct(pU,pB_pU);
        double c=vectorNorm2(pU)-delta*delta;
        double k1=(-b-sqrt(b*b-4*a*c))/(2*a)+1;
        double k2=(-b+sqrt(b*b-4*a*c))/(2*a)+1;
        double k=fmax(k1,k2);

        if (0<=k && k<=1)
            vectorScalarProduct(pC,pU,k);

        else if (1<=k && k<=2){
            vectorScalarProduct(pC,pB_pU,k-1);
            vectorAdd(pU,pC,pC);
        }
        FreeVectorMem(pB_pU);
        cout<<"tao ";
    }

    FreeVectorMem(pU);
    FreeVectorMem(pB);
    FreeVectorMem(gTB);
    FreeVectorMem(aux);
}


#endif // DOGLEG_H_INCLUDED
