#ifndef ALPHABACK_H_INCLUDED
#define ALPHABACK_H_INCLUDED

#include "MatrixInterface.h"

double InertiaBacktracking(condition fc,double alpha1,Vector *x, funcion fx, Vector*p, Vector*grad,Matrix *Hess)
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
    double alpha2=(alpha1*alpha1*phi02)/(2*(phi01+phi02*alpha1-phi_alpha1));

    int maxite=10,count=0,tracking=0,x1=2,x2=5,x3=3;
    double alphaB=alpha2;
    while(fc(alphaB,x,fx,p,grad)==0 && count<maxite)
    {
        tracking++;
        if (tracking>x2)
        {
            alphaB=x3*alphaB;
            tracking=0;
        }
        else
            alphaB=alphaB/x1;
        count++;
    }

    FreeVectorMem(alpha_p);
    FreeVectorMem(x_new);

    if (count==10)
        return alpha2;

    else
        return alphaB;
}

#endif // ALPHABACK_H_INCLUDED
