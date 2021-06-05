#ifndef GRAD_H_INCLUDED
#define GRAD_H_INCLUDED

#include "MatrixInterface.h"

extern double c1,c2;
void GradientDescent(Vector *x,funcion fx,Vector *grad, Matrix *Hess,Vector *p)
{
    c1=1e-4;c2=0.1;
    Gradient(fx,x,p);
    for (int i=0;i<dim;i++)
        p->data[i]=-p->data[i];
}

#endif // GRAD_H_INCLUDED
