#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "MatrixInterface.h"
#include "PGMFunctions.h"

using namespace std;

///Funcion que hace la rotacion a las coordenadas
///xRot contiene las coordenadas de la imagen original con las cuales interpolaremos
void CoordTransformation(Vector *x,Vector *alpha, Vector *xRot)
{
    xRot->data[0]=x->data[0]+x->data[0]*alpha->data[0]+x->data[1]*alpha->data[1]+alpha->data[2];
    xRot->data[1]=x->data[1]+x->data[0]*alpha->data[3]+x->data[1]*alpha->data[4]+alpha->data[5];
}

///funcion que interpola el valor de la imagen en el pixel rotado
double BilinearInterpolation(Vector *xRot, Matrix *Image)
{
    double delta1=xRot->data[0]-floor(xRot->data[0]);
    double delta2=xRot->data[1]-floor(xRot->data[1]);

    int j=(int)xRot->data[0]-imgsize/2;
    int i=imgsize/2-(int)xRot->data[1];

    double interpolation=0;

    interpolation=(1-delta1)*(1-delta2)*Image->data[((i+2*imgsize)%imgsize*imgsize)+(j+2*imgsize)%imgsize]
                        +(1-delta1)*delta2*Image->data[((i+2*imgsize)%imgsize*imgsize)+((j+1)+2*imgsize)%imgsize]
                        +delta1*(1-delta2)*Image->data[(((i+1)+2*imgsize)%imgsize*imgsize)+(j+2*imgsize)%imgsize]
                        +delta1*delta2*Image->data[(((i+1)+2*imgsize)%imgsize*imgsize)+((j+1)+2*imgsize)%imgsize];

    return floor(interpolation);

}

///calcula el producto gradI1(xRot)*X
void GradXProduct(Matrix **ImgGrad,Vector *x,Vector *xRot, Vector *Product)
{
    Vector *grad=RequestVectorDoubleMem(2);

    grad->data[0]=BilinearInterpolation(xRot,ImgGrad[0]);
    grad->data[1]=BilinearInterpolation(xRot,ImgGrad[1]);

    Product->data[0]=x->data[0]*grad->data[0];
    Product->data[1]=x->data[1]*grad->data[0];
    Product->data[2]=grad->data[0];

    Product->data[3]=x->data[0]*grad->data[1];
    Product->data[4]=x->data[1]*grad->data[1];
    Product->data[5]=grad->data[1];

    FreeVectorMem(grad);
}

double GradXalphaProduct(Matrix **ImgGrad, Vector *x, Vector *xRot, Vector *alpha)
{
    Vector *grad=RequestVectorDoubleMem(2);
    Vector *xalpha=RequestVectorDoubleMem(2);

    grad->data[0]=BilinearInterpolation(xRot,ImgGrad[0]);
    grad->data[1]=BilinearInterpolation(xRot,ImgGrad[1]);

    xalpha->data[0]=x->data[0]*alpha->data[0]+x->data[1]*alpha->data[1]+alpha->data[2];
    xalpha->data[1]=x->data[0]*alpha->data[3]+x->data[1]*alpha->data[4]+alpha->data[5];

    double product=vectorInnerProduct(grad,xalpha);

    FreeVectorMem(grad);
    FreeVectorMem(xalpha);

    return product;
}

void PrintImage(Matrix *Image1,Vector *alpha, Matrix *Image2,char *filename)
{
    Vector *x=RequestVectorDoubleMem(2);
    Vector *xRot=RequestVectorDoubleMem(2);

    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
        {
            ///guardamos en el vector x las coordenadas que vamos a utilizar
                x->data[0]=j-imgsize/2;
                x->data[1]=imgsize/2-i;

            ///transformacion inversa de coordenadas de acuerdo al vector alpha
            CoordTransformation(x,alpha,xRot);

            Image2->data[i*imgsize+j]=BilinearInterpolation(xRot,Image1);
        }
    writePGM(filename,Image2);
    FreeVectorMem(x);
    FreeVectorMem(xRot);
}

///funcion que calcula el error que la rotacion alpha produce en la imagen observada con respecto a la de referencia
///recibe la imagen original, sin rotar, su gradiente, y la imagen de referencia
double ErrorFunction(Matrix *Image1, Matrix **ImgGrad, Matrix *ImageRef, Vector *alpha, Vector *pk)
{
    double error=0;
    Vector *x=RequestVectorDoubleMem(2);
    Vector *xRot=RequestVectorDoubleMem(2);

    ///barremos la imagen pixel por pixel
    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
        {
            ///guardamos en el vector x las coordenadas que vamos a utilizar
                x->data[0]=j-imgsize/2;
                x->data[1]=imgsize/2-i;

            ///transformacion de coordenadas de acuerdo al vector alpha
            CoordTransformation(x,alpha,xRot);

            ///I2(x)-I1(xRot)
            double sum=ImageRef->data[i*imgsize+j];
            sum-=BilinearInterpolation(xRot,Image1);
            sum-=GradXalphaProduct(ImgGrad,x,xRot,pk);

            ///se acumula el error por pixel
            error+=0.5*sum*sum;
        }
    FreeVectorMem(x);
    FreeVectorMem(xRot);
    return error;
}

///Calcula el gradiente de la funcion de error, un vector de tamanio "dim"
void ErrorGradient(Matrix *Image1, Matrix **ImgGrad, Matrix *ImageRef,Vector *alpha, Vector *gradient, Vector *pk)
{
    vectorFill(gradient,0);
    Vector *x=RequestVectorDoubleMem(2);
    Vector *xRot=RequestVectorDoubleMem(2);
    Vector *Product=RequestVectorDoubleMem(dim);

        ///barremos la imagen pixel por pixel
        for (int i=0;i<imgsize;i++)
            for (int j=0;j<imgsize;j++)
            {
                ///guardamos en el vector x las coordenadas que vamos a utilizar
                x->data[0]=j-imgsize/2;
                x->data[1]=imgsize/2-i;

                ///transformacion de coordenadas de acuerdo al vector alpha
                CoordTransformation(x,alpha,xRot);

                ///se calcula el escalar -I2(x)+I1(xRot)
                double scalar=BilinearInterpolation(xRot,Image1);
                scalar-=ImageRef->data[i*imgsize+j];
                scalar+=GradXalphaProduct(ImgGrad,x,xRot,pk);

                ///calcula el producto gradI1(xRot)*X
                GradXProduct(ImgGrad,x,xRot,Product);
                vectorScalarProduct(Product,Product,scalar);

                ///se acumula el resultado
                vectorAdd(Product,gradient,gradient);
            }
    FreeVectorMem(x);
    FreeVectorMem(xRot);
    FreeVectorMem(Product);
}

///Calcula el hessiano de la funcion de error, una matrix de dimension "dim"
void ErrorHessian(Matrix *Image1, Matrix **ImgGrad, Matrix **ImgHessian, Matrix *ImageRef,Vector *alpha, Matrix *Hessian)
{
    matrixFill(Hessian,0);
    Vector *x=RequestVectorDoubleMem(2);
    Vector *xRot=RequestVectorDoubleMem(2);
    Vector *Product1=RequestVectorDoubleMem(dim);
    Matrix *Product2=RequestMatrixDoubleMem(dim,dim);
    Matrix *Product3=RequestMatrixDoubleMem(dim,dim);

        ///barremos la imagen pixel por pixel
        for (int i=0;i<imgsize;i++)
            for (int j=0;j<imgsize;j++)
            {
                ///guardamos en el vector x las coordenadas que vamos a utilizar
                x->data[0]=j-imgsize/2;
                x->data[1]=imgsize/2-i;

                ///transformacion de coordenadas de acuerdo al vector alpha
                CoordTransformation(x,alpha,xRot);

                ///calcula el producto gradI1(xRot)*X
                GradXProduct(ImgGrad,x,xRot,Product1);

                ///calcula el producto gradI1(xRot)*X(gradI1(xRot)*X)^T
                vectorOuterProduct(Product1,Product1,Product2);

                ///se acumula el resultado
                matrixAdd(Product2,Hessian,Hessian);
            }
    FreeVectorMem(x);
    FreeVectorMem(xRot);
    FreeVectorMem(Product1);
    FreeMatrixMem(Product2);
    FreeMatrixMem(Product3);
}

#endif // FUNCTIONS_H_INCLUDED
