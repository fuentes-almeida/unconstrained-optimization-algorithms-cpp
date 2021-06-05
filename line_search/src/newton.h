#ifndef NEWTON_H_INCLUDED
#define NEWTON_H_INCLUDED

#include "MatrixInterface.h"

extern double c1,c2;
void NewtonMethod(Vector *x,funcion fx,Vector *grad, Matrix *Hess,Vector *p)
{
    c1=1e-4;c2=0.9;
    Matrix *HessCopy=RequestMatrixDoubleMem(dim,dim);
    for (int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            HessCopy->data[i*dim+j]=Hess->data[i*dim+j];

    choleskySolve(HessCopy,grad,p);
        for (int i=0;i<dim;i++)
        p->data[i]=-p->data[i];

    FreeMatrixMem(HessCopy);
}

#endif // NEWTON_H_INCLUDED


/*


*/
