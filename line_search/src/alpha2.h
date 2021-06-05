#ifndef ALPHA2_H_INCLUDED
#define ALPHA2_H_INCLUDED

#include "MatrixInterface.h"

double CuadraticModel(condition fc,double alpha1,Vector *x, funcion fx, Vector*p, Vector*grad,Matrix *Hess)
{
    if (fc(alpha1,x,fx,p,grad)==1)
        return alpha1;

    double phi01=fx(x);
    double phi02=vectorInnerProduct(grad,p);

    Vector *alpha_p=RequestVectorDoubleMem(dim);
    Vector *x_new=RequestVectorDoubleMem(dim);

    vectorScalarProduct(alpha_p,p,alpha1);
    vectorAdd(x,alpha_p,x_new);

    double phi_alpha1=fx(x_new);
    double alpha2=(alpha1*alpha1*phi02)/(2.0*(phi01+phi02*alpha1-phi_alpha1));

    FreeVectorMem(alpha_p);
    FreeVectorMem(x_new);

    return alpha2;
}

#endif // ALPHA2_H_INCLUDED
