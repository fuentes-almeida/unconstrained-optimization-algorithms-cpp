#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Configuration.h"
#include "MatrixInterface.h"
#include "newton.h"
#include "grad.h"
#include "wolfeSmooth.h"
#include "wolfeStrong.h"
#include "goldstein.h"
#include "alpha2.h"
#include "alpha3.h"
#include "alphaback.h"

void MinimizeFunction(Vector *x,double alphak,funcion fx,method fm,condition fc,model fa)
{
    clock_t start=clock();
    int max_ite=1e6; double tole=1e-6;
    Vector *grad=RequestVectorDoubleMem(dim);
    Vector *p=RequestVectorDoubleMem(dim);
    Vector *xnew=RequestVectorDoubleMem(dim);
    Matrix *Hess=RequestMatrixDoubleMem(dim,dim);
    double alpha1;

    int count=0;
    double dist=1000;
    while (count<max_ite && dist>tole)
    {
        Gradient(fx,x,grad);
        Hessian(fx,x,Hess);

        fm(x,fx,grad,Hess,p);
        alpha1=fa(fc,alphak,x,fx,p,grad,Hess);

        for (int i=0;i<dim;i++)
            xnew->data[i]=x->data[i]+alpha1*p->data[i];
        alphak=alpha1;

        dist=Euclidean(x,xnew);
            for (int i=0;i<dim;i++)
                x->data[i]=xnew->data[i];
        count++;
    }
    clock_t end=clock();

    printf("%d %lf %lf\n",count, 1000*(end-start)/(double)CLOCKS_PER_SEC,fx(x));

    FreeMatrixMem(Hess);
    FreeVectorMem(grad);
    FreeVectorMem(p);
    FreeVectorMem(xnew);
}

int main(int argc, char** argv)
{
    srand(time(NULL));
    if (argc!=4){printf("Missing Arguments");exit(-1);}

    int method_ref=atoi(argv[1]);
    int cond=atoi(argv[2]);
    int alpha_mode=atoi(argv[3]);

    method *M=(method*)malloc(2*sizeof(method));
    M[0]=GradientDescent;
    M[1]=NewtonMethod;

    condition *C=(condition*)malloc(3*sizeof(condition));
    C[0]=SmoothWolfe;
    C[1]=StrongWolfe;
    C[2]=Goldstein;

    model *A=(model*)malloc(3*sizeof(model));
    A[0]=CuadraticModel;
    A[1]=CubicModel;
    A[2]=InertiaBacktracking;

    double alpha0=0.2;
    Vector *x=RequestVectorDoubleMem(dim);

    for (int k=0;k<30;k++)
    {
        for (int i=0;i<dim;i++)
            x->data[i]=15.0*rand()/RAND_MAX-5.0;
        MinimizeFunction(x,alpha0,Rosenbrock,M[method_ref],C[cond],A[alpha_mode]);
    }

    free(M);
    free(C);
    free(A);
    FreeVectorMem(x);

    return 0;
}
