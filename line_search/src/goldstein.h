#ifndef GOLDSTEIN_H_INCLUDED
#define GOLDSTEIN_H_INCLUDED

#include "MatrixInterface.h"

int Goldstein(double alpha,Vector *x,funcion fx,Vector *p,Vector *grad)
{
    double c=0.05;
    Vector *alpha_p=RequestVectorDoubleMem(dim);
    Vector *x_new=RequestVectorDoubleMem(dim);

    vectorScalarProduct(alpha_p,p,alpha);
    vectorAdd(x,alpha_p,x_new);
    double fx2=fx(x_new);


    double grad_x_p=vectorInnerProduct(grad,p);
    double aux1=fx(x)+c*alpha*grad_x_p;
    double aux2=fx(x)+(1-c)*alpha*grad_x_p;

    FreeVectorMem(alpha_p);
    FreeVectorMem(x_new);

    return (aux2<=fx2 && fx2<=aux1);
}

#endif // GOLDSTEIN_H_INCLUDED
