#ifndef WOLFESTRONG_H_INCLUDED
#define WOLFESTRONG_H_INCLUDED

#include "MatrixInterface.h"
extern double c1,c2;
int StrongWolfe(double alpha,Vector *x,funcion fx,Vector *p,Vector *grad)
{
    Vector *alpha_p=RequestVectorDoubleMem(dim);
    Vector *x_new=RequestVectorDoubleMem(dim);
    Vector *grad_xnew=RequestVectorDoubleMem(dim);

    vectorScalarProduct(alpha_p,p,alpha);
    vectorAdd(x,alpha_p,x_new);


    double grad_x_p=vectorInnerProduct(grad,p);
    double aux1=fx(x)+c1*alpha*grad_x_p;

    Gradient(fx,x_new,grad_xnew);

    double aux2=vectorInnerProduct(grad_xnew,p);

    FreeVectorMem(alpha_p);
    FreeVectorMem(x_new);
    FreeVectorMem(grad_xnew);

    return (fx(x_new)<=aux1 && fabs(aux2)<=c2*fabs(grad_x_p));
}

#endif // WOLFESTRONG_H_INCLUDED
