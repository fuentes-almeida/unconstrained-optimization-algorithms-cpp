#ifndef MRC_H_INCLUDED
#define MRC_H_INCLUDED

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "MatrixInterface.h"
#include "PGMFunctions.h"
#include "dogleg.h"
#include "functions.h"
using namespace std;

///funcion que aproxima el error por medio de la serie de taylor
///mk(p) = fk + grad(fk)^T*p + 0.5p^T*Bk*p
double Model(Vector *p, Vector *gradient, Matrix *Hessian)
{
    Vector *aux=RequestVectorDoubleMem(dim);
    ///segundo termino
    double Gradfk_p=vectorInnerProduct(p,gradient);

    ///tercer termino
    matrixVectorProduct(Hessian,p,aux);
    double pTBp=vectorInnerProduct(aux,p);
    FreeVectorMem(aux);

    return -Gradfk_p-0.5*pTBp;
}

///Calcular la region de confianza: deltaMax > 0, 0 < delta0 < deltaMax, 0 <= threshold < 0.25
void TrustRegion(Matrix *Image1, Matrix **ImgGrad, Matrix **ImgHessian, Matrix *ImageRef)
{
    double deltaMax=1, delta=0.5, epsilon=1000, threshold=0.1, fk_prev=0, fk_new=0; int k=1;

    Matrix *ImageAux=RequestMatrixDoubleMem(imgsize,imgsize);
    Vector *gradient=RequestVectorDoubleMem(dim);
    Matrix *Hessian=RequestMatrixDoubleMem(dim,dim);
    Vector *pk=RequestVectorDoubleMem(dim);
    Vector *alphaNew=RequestVectorDoubleMem(dim);
    Vector *aux=RequestVectorDoubleMem(dim);

    Vector *alpha=RequestVectorDoubleMem(dim);
    double theta=0;
    alpha->data[0]=cos(theta)-1;
    alpha->data[1]=sin(theta);
    alpha->data[2]=10;
    alpha->data[3]=-alpha->data[1];
    alpha->data[4]=alpha->data[0];
    alpha->data[5]=-10;

    string filename="output0.pgm";
    PrintImage(Image1,alpha,ImageAux,&filename[0]);

    fk_prev=ErrorFunction(Image1,ImgGrad,ImageRef,alpha,pk);

    while (fk_prev>epsilon && k<300)
    {
        ErrorGradient(Image1,ImgGrad,ImageRef,alpha,gradient,pk);
        ErrorHessian(Image1,ImgGrad,ImgHessian,ImageRef,alpha,Hessian);

        Dogleg(delta,gradient,Hessian,pk);

        fk_new=ErrorFunction(Image1,ImgGrad,ImageRef,alpha,pk);
        double mk=Model(pk,gradient,Hessian);
        double radius=(fk_prev-fk_new)/mk;

        if (radius<0.25)
            delta=0.25*delta;
        else if (radius>0.75 && fabs(vectorNorm(pk)-delta)<0.01)
            delta=fmin(2*delta,deltaMax);

        if (radius>threshold)
        {
            vectorAdd(alpha,pk,alpha);
            fk_prev=fk_new;

            string outputfile="output"+to_string(k)+".pgm";
            PrintImage(Image1,alpha,ImageAux,&outputfile[0]);
            k++;
        }

        double product=vectorInnerProduct(gradient,pk);
        cout<<setprecision(10)<<product<<" "<<radius<<" "<<fk_prev<<" "<<endl;
    }

    FreeMatrixMem(ImageAux);
    FreeMatrixMem(Hessian);
    FreeVectorMem(gradient);
    FreeVectorMem(pk);
    FreeVectorMem(alphaNew);
    FreeVectorMem(aux);
}

#endif // MRC_H_INCLUDED
