#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "CommonTypes.h"

/*!
\brief Modulo exclusivo para contener variables
globales de configuracion que se emplean en los
diferentes metodos de optimizacion y que se inicializan
dinamicamente en cada problema a resolver.

Cada variable al agregarse debe tener la siguiente notacion g_nombreVariable
seguida de un comentario breve sobre el uso de la variable.
*/
double c1,c2;

struct Config
{
	int maxIter; //Maximo de iteraciones a ejecutar en cada metodo
	double eps; //Epsilon del sistema
};

extern Config g_config;

#endif

/*
///calcula el producto gradI1(xRot)*X*alpha
double GradXalphaProduct(Matrix **ImgGrad, Vector *x, Vector *xRot, Vector *alpha)
{
    Vector *grad=RequestVectorDoubleMem(2);
    Vector *xalpha=RequestVectorDoubleMem(2);

    grad->data[0]=BilinearInterpolation(xRot,ImgGrad[0]);
    grad->data[1]=BilinearInterpolation(xRot,ImgGrad[1]);

    //xalpha->data[0]=alpha->data[0]-alpha->data[1]*x->data[1];
   // xalpha->data[1]=alpha->data[2]+alpha->data[1]*x->data[0];

    xalpha->data[0]=x->data[0]*alpha->data[0]+x->data[1]*alpha->data[1]+alpha->data[2];
    xalpha->data[1]=x->data[0]*alpha->data[3]+x->data[1]*alpha->data[4]+alpha->data[5];

    double product=vectorInnerProduct(grad,xalpha);

    FreeVectorMem(grad);
    FreeVectorMem(xalpha);

    return product;
}


void BFGS(Vector *alpha,Vector *alphaNew, Vector *Gradient, Vector *GradientNew,Matrix *Hessian)
{
    Vector *Yk=RequestVectorDoubleMem(dim);
    Vector *Sk=RequestVectorDoubleMem(dim);
    Matrix *SkSkT=RequestMatrixDoubleMem(dim,dim);
    Matrix *Product1=RequestMatrixDoubleMem(dim,dim);
    Matrix *Product2=RequestMatrixDoubleMem(dim,dim);
    Matrix *YkYkT=RequestMatrixDoubleMem(dim,dim);

    vectorSubstract(GradientNew,Gradient,Yk);
    vectorSubstract(alphaNew,alpha,Sk);

    vectorOuterProduct(Sk,Sk,SkSkT);
    matrixProduct(Hessian,SkSkT,Product1);
    matrixProduct(Product1,Hessian,Product2);
    double SkBkSk=quadraticFormProduct(Hessian,Sk);
    matrixScalarProduct(Product2,Product2,-1/SkBkSk);

    vectorOuterProduct(Yk,Yk,YkYkT);
    double YkSk=vectorInnerProduct(Yk,Sk);
    matrixScalarProduct(YkYkT,YkYkT,1/YkSk);

    matrixSubstract(Hessian,Product2,Hessian);
    matrixAdd(Hessian,YkYkT,Hessian);

}

double Model2(Matrix *Image1, Matrix **ImgGrad, Matrix *ImageRef, Vector *alpha)
{
    double error=0;
    Vector *x=RequestVectorDoubleMem(2);
    Vector *xRot=RequestVectorDoubleMem(2);

    ///barremos la imagen pixel por pixel
    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
        {
            ///guardamos en el vector x las coordenadas que vamos a utilizar
            x->data[0]=j;
            x->data[1]=i;

            ///transformacion de coordenadas de acuerdo al vector alpha
            CoordTransformation(x,alpha,xRot);

            ///I2(x)-I1(xRot)
            double sum=ImageRef->data[i*imgsize+j];
            sum-=BilinearInterpolation(xRot,Image1);
            sum-=GradXalphaProduct(ImgGrad,x,xRot,alpha);

            ///se acumula el error por pixel
            error+=0.5*sum*sum;
        }
    FreeVectorMem(x);
    FreeVectorMem(xRot);
    return error;
}

void HessianProduct(Vector *x, Vector *xRot, Matrix **ImgHessian, Matrix *Product)
{
    Matrix *hess=RequestMatrixDoubleMem(2,2);
    for (int i=0;i<4;i++)
        hess->data[i]=BilinearInterpolation(xRot,ImgHessian[i]);

    Matrix *X=RequestMatrixDoubleMem(dim,2);
    Matrix *XT=RequestMatrixDoubleMem(2,dim);
        X->data[0]=X->data[7]=x->data[0];
        X->data[2]=X->data[9]=x->data[1];
        X->data[4]=X->data[11]=1;
    matrixTranspose(X,XT);

    Matrix *XThess=RequestMatrixDoubleMem(dim,2);
    matrixProduct(XT,hess,XThess);
    matrixProduct(XThess,X,Product);

    FreeMatrixMem(X);
    FreeMatrixMem(XT);
    FreeMatrixMem(XThess);
    FreeMatrixMem(hess);
}



double Model2(Matrix *Image1, Matrix **ImgGrad, Matrix *ImageRef, Vector *alpha)
{
    double error=0;
    Vector *x=RequestVectorDoubleMem(2);
    Vector *xRot=RequestVectorDoubleMem(2);

    ///barremos la imagen pixel por pixel
    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
        {
            ///guardamos en el vector x las coordenadas que vamos a utilizar
            x->data[0]=j;
            x->data[1]=i;

            ///transformacion de coordenadas de acuerdo al vector alpha
            CoordTransformation(x,alpha,xRot);

            ///I2(x)-I1(xRot)
            double sum=ImageRef->data[i*imgsize+j];
            sum-=BilinearInterpolation(xRot,Image1);
            sum-=GradXalphaProduct(ImgGrad,x,xRot,alpha);

            ///se acumula el error por pixel
            error+=0.5*sum*sum;
        }
    FreeVectorMem(x);
    FreeVectorMem(xRot);
    return error;
}

void UpdateAlpha(Vector *alpha, Vector *pk, Vector *alphaNew)
{
    alphaNew->data[0]=alpha->data[0]+pk->data[0];
    alphaNew->data[1]=alpha->data[1]+pk->data[1];
    alphaNew->data[3]=-alphaNew->data[1];
    alphaNew->data[4]=alpha->data[0]+pk->data[0];
    alphaNew->data[2]=alpha->data[2]+pk->data[2];
    alphaNew->data[5]=alpha->data[5]+pk->data[5];
}
*/
