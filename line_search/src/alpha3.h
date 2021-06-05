#ifndef ALPHA3_H_INCLUDED
#define ALPHA3_H_INCLUDED

#include "MatrixInterface.h"

double CubicModel(condition fc,double alpha1,Vector *x, funcion fx, Vector*p, Vector*grad,Matrix *Hess)
{
    if (fc(alpha1,x,fx,p,grad)==1)
        return alpha1;

    Vector *product=RequestVectorDoubleMem(dim);
    Vector *alpha_p=RequestVectorDoubleMem(dim);
    Vector *x_new=RequestVectorDoubleMem(dim);

    matrixVectorProduct(Hess,p,product);

    double phi01=fx(x);
    double phi02=vectorInnerProduct(grad,p);
    double phi03=vectorInnerProduct(p,product);

    double b=phi02;
    double c=0.5*phi03;

    vectorScalarProduct(alpha_p,p,alpha1);
    vectorAdd(x,alpha_p,x_new);

    double phi_alpha2=fx(x_new);
    double d=(phi_alpha2-phi01-alpha1*phi02-0.5*alpha1*alpha1*phi03)/(alpha1*alpha1*alpha1);

    if (c*c-3*b*d<0) return alpha1;

    double alpha3=(-c+sqrt(c*c-3*b*d))/(3*d);

    FreeVectorMem(alpha_p);
    FreeVectorMem(x_new);
    FreeVectorMem(product);

    if (alpha3>0) return alpha3;
    else return alpha1;

}

#endif // ALPHA3_H_INCLUDED
