#ifndef GCNL_H_INCLUDED
#define GCNL_H_INCLUDED

#include "function.h"
#include "lineSearch.h"

///calcula la beta de acuerdo al indicador dado como parametro
double CalculateBeta(Matrix *gradk, Matrix *gradk_new, Matrix *pk, int flag)
{
    double beta=0,gkgk,gkgk_new;
    int N=gradk->rows, M=gradk->cols;
    Matrix *Aux=RequestMatrixDoubleMem(N,M);

    switch (flag)
    {
        case 1:///beta de Fletcher-Reeves
        gkgk=matrixInnerProduct(gradk,gradk);
        gkgk_new=matrixInnerProduct(gradk_new,gradk_new);
        beta=gkgk_new/gkgk; break;

        case 2:///beta de Polak-Ribiere
        gkgk=matrixInnerProduct(gradk,gradk);
        matrixSubstract(gradk_new,gradk,Aux);
        gkgk_new=matrixInnerProduct(gradk_new,Aux);
        beta=gkgk_new/gkgk;
        if (beta<0) beta=0; break;

        case 3:///beta de Hestenes-Steifel
        matrixSubstract(gradk_new,gradk,Aux);
        gkgk_new=matrixInnerProduct(gradk_new,Aux);
        gkgk=matrixInnerProduct(pk,Aux);
        beta=gkgk_new/gkgk; break;
    }

    FreeMatrixMem(Aux);
    return beta;
}

///Suavizamiento de la imagen utilizando el metodo de Gradiente Conjugado No Lineal (con HQR)
void NonLinearConjugatedGradient(Matrix* ImageX, Matrix* ImageY, double lambda, double k_Const, int MD)
{
    int k=0, Ite_max=500, N=ImageX->rows, M=ImageX->cols;
    double epsilon=100, alpha=0.1, error=1000;
    const char *outputfile="NLConjGrad.pgm";

    Matrix *gradk=RequestMatrixDoubleMem(N,M);
    Matrix *gradk_new=RequestMatrixDoubleMem(N,M);
    Matrix *pk=RequestMatrixDoubleMem(N,M);
    Matrix *alpha_pk=RequestMatrixDoubleMem(N,M);

    ///valores iniciales del algoritmo asumiendo ws = 1 para todos los vecinos
    double fx=ModelFunction(ImageX,ImageY,lambda,k_Const,1);       ///f0 = f(x0)
    ModelFunctionGradient(ImageX,ImageY,gradk,lambda,k_Const,1);   ///g0 = grad(f0)
    matrixScalarProduct(pk,gradk,-1);                              ///p0 = -g0

    while (fabs(error)>epsilon && k<Ite_max)   ///condicion de paro fijada por el modulo del gradiente
    {
        error=fx; k++;
        ///alpha se calcula por busqueda en linea, con un modelo cuadratico y condicion fuerte de Wolfe
        alpha=CuadraticModel(alpha,ImageX,ImageY,pk,gradk,fx,lambda,k_Const);
         printf("%d %lf %lf ",k, fx, alpha);

        ///calculo de la nueva imagen en funcion de alpha y pk
        matrixScalarProduct(alpha_pk,pk,alpha);
        matrixAdd(ImageX,alpha_pk,ImageX);

        ///calculo del nuevo gradiente, tomando en cuenta ws en funcion de la nueva imagen
        ModelFunctionGradient(ImageX,ImageY,gradk_new,lambda,k_Const,0);

        ///calculo de beta
        double beta=CalculateBeta(gradk,gradk_new,pk,MD);
        printf("%lf ",beta);

        ///calculo del nuevo paso pk+1
        matrixScalarProduct(pk,pk,beta);
        matrixSubstract(pk,gradk_new,pk);

        ///calculo del nuevo valor de la funcion objetivo
        fx=ModelFunction(ImageX,ImageY,lambda,k_Const,0);

        ///se sustituye el vejo gradiente por el nuevo
        Matrix_copy(gradk_new,gradk);
        printf("%lf\n",matrixNorm(gradk));

        ///calculo del nuevo error entre imagenes
        error=error-fx;
    }

    NormalizeImage(ImageX,ImageX);
    writePGM(outputfile,ImageX);

    FreeMatrixMem(gradk);
    FreeMatrixMem(gradk_new);
    FreeMatrixMem(pk);
    FreeMatrixMem(alpha_pk);
}

#endif // GCNL_H_INCLUDED
