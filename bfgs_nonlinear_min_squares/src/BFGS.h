#ifndef BFGS_H_INCLUDED
#define BFGS_H_INCLUDED

#include "functions.h"

using namespace std;

void BFGS_step(double *fk_cur, Vector *Rx, Vector *Xk_cur, Vector *Gk_cur, Matrix *Hk_cur, Matrix *Jacobian, Vector *Alpha, Vector *Mu, Vector *Sigma, Matrix *ImgWrapped, Matrix **ImgWrappedGrad, int size, int flag, int methodID, int *index)
{
    Vector *Pk=RequestVectorDoubleMem(size);
    Vector *Xk_new=RequestVectorDoubleMem(size);
    Vector *Gk_new=RequestVectorDoubleMem(size);
    Matrix *Hk_new=RequestMatrixDoubleMem(size,size);

    Vector *Yk=RequestVectorDoubleMem(size);
    Vector *Sk=RequestVectorDoubleMem(size);

    Matrix *Aux1=RequestMatrixDoubleMem(size,size);
    Matrix *Aux2=RequestMatrixDoubleMem(size,size);
    Matrix *Aux3=RequestMatrixDoubleMem(size,size);

        ///Calculo de la direccion de descenso Pk
        matrixVectorProduct(Hk_cur,Gk_cur,Pk);
        vectorScalarProduct(Pk,Pk,-1);

        ///calculo del tamanio de paso alfa
        double StepSize=LineSearch(1,*fk_cur,Xk_cur,Gk_cur,Pk,ImgWrapped,ImgWrappedGrad,Alpha,Mu,Sigma,size,flag,methodID,index);

        ///calculo del nuevo vector de X
        vectorScalarProduct(Pk,Pk,StepSize);
        vectorAdd(Xk_cur,Pk,Xk_new);

        if (flag==1)
            vectorCopy(Xk_new,Alpha);
        if (flag==2)
            vectorCopy(Xk_new,Mu);
        if (flag==3)
            vectorCopy(Xk_new,Sigma);

        ///Calculo del nuevo gradiente
        double fk_new=ObjectiveFunction(Rx,ImgWrappedGrad,Alpha,Mu,Sigma);
        ConstructJacobian(Jacobian,Alpha,Mu,Sigma,flag,methodID,index);
        matrixVectorProductRandom(Jacobian,Rx,Gk_new,methodID,index);

        ///calculo de Sk y Yk
        vectorSubstract(Xk_new,Xk_cur,Sk);
        vectorSubstract(Gk_new,Gk_cur,Yk);

        ///calculo de RHo
        double Rho=1/vectorInnerProduct(Sk,Yk);

        ///Calculo de la nuevo hessiano inverso
        vectorOuterProduct(Sk,Yk,Aux1);
        matrixScalarProduct(Aux1,Aux1,-Rho);
        matrixAddToDiagonal(Aux1,Aux1,1);

        vectorOuterProduct(Yk,Sk,Aux2);
        matrixScalarProduct(Aux2,Aux2,-Rho);
        matrixAddToDiagonal(Aux2,Aux2,1);

        matrixProduct(Aux1,Hk_cur,Aux3);
        matrixProduct(Aux3,Aux2,Hk_new);

        vectorOuterProduct(Sk,Sk,Aux3);
        matrixScalarProduct(Aux3,Aux3,Rho);
        matrixAdd(Hk_new,Aux3,Hk_new);

        ///Actualizacion de variables
        *fk_cur=fk_new;
        vectorCopy(Xk_new,Xk_cur);
        vectorCopy(Gk_new,Gk_cur);
        Matrix_copy(Hk_new,Hk_cur);

    FreeVectorMem(Pk);
    FreeVectorMem(Sk);
    FreeVectorMem(Yk);
    FreeMatrixMem(Aux1);
    FreeMatrixMem(Aux2);
    FreeMatrixMem(Aux3);
    FreeMatrixMem(Hk_new);
    FreeVectorMem(Xk_new);
    FreeVectorMem(Gk_new);
}

void NonLinearLeastSquare(Matrix *ImgWrapped, int MethodID)
{
    int N=imgsize, M=imgsize, k=0, size, *index;
    double  fk_cur, fk_new, Norm=1000;

    if (MethodID==1) size=N*M;
    else { size=N*M/2; index=(int*)calloc(size,sizeof(int)); }

    Matrix **ImgWrappedGrad=(Matrix**)malloc(2*sizeof(Matrix*));
    CalculateWrappingGradients(ImgWrapped,ImgWrappedGrad);

    Matrix *ImageX=RequestMatrixDoubleMem(N,M);
    Vector *Rx=RequestVectorDoubleMem(2*N*M);

    Vector *Alpha=RequestVectorDoubleMem(Nnum); vectorFill(Alpha,-1);
    Vector *Mu=RequestVectorDoubleMem(2*Nnum);  CalculateGaussianMeans(Mu,N,M);
    Vector *Sigma=RequestVectorDoubleMem(Nnum); vectorFill(Sigma,imgsize/(Num));

    Matrix *Jacobian1=RequestMatrixDoubleMem(Nnum,2*size);  Matrix *Jacobian2=RequestMatrixDoubleMem(2*Nnum,2*size);    Matrix *Jacobian3=RequestMatrixDoubleMem(Nnum,2*size);
    Matrix *Hk_cur1=RequestMatrixDoubleMem(Nnum,Nnum);      Matrix *Hk_cur2=RequestMatrixDoubleMem(2*Nnum,2*Nnum);      Matrix *Hk_cur3=RequestMatrixDoubleMem(Nnum,Nnum);
    Vector *Gk_cur1=RequestVectorDoubleMem(Nnum);           Vector *Gk_cur2=RequestVectorDoubleMem(2*Nnum);             Vector *Gk_cur3=RequestVectorDoubleMem(Nnum);
    Vector *Xk_cur1=RequestVectorDoubleMem(Nnum);           Vector *Xk_cur2=RequestVectorDoubleMem(2*Nnum);             Vector *Xk_cur3=RequestVectorDoubleMem(Nnum);

    UnWrapImage(Alpha,Mu,Sigma,ImageX);
    NormalizeImage(ImageX,ImageX,255);
    writePGM("Result0.pgm",ImageX);

    while (k<50)
    {
        if (MethodID==2) RandomVector(size,N*M,index);

        if (k==0)
        {
            ///Valores iniciales
            matrixAddToDiagonal(Hk_cur1,Hk_cur1,0.5);
            matrixAddToDiagonal(Hk_cur2,Hk_cur2,0.5);
            matrixAddToDiagonal(Hk_cur3,Hk_cur3,0.5);

            vectorCopy(Alpha,Xk_cur1);
            vectorCopy(Mu,Xk_cur2);
            vectorCopy(Sigma,Xk_cur3);

            fk_new=ObjectiveFunction(Rx,ImgWrappedGrad,Alpha,Mu,Sigma);

            ConstructJacobian(Jacobian1,Alpha,Mu,Sigma,1,MethodID,index);
            ConstructJacobian(Jacobian2,Alpha,Mu,Sigma,2,MethodID,index);
            ConstructJacobian(Jacobian3,Alpha,Mu,Sigma,3,MethodID,index);

            matrixVectorProductRandom(Jacobian1,Rx,Gk_cur1,MethodID,index);
            matrixVectorProductRandom(Jacobian2,Rx,Gk_cur2,MethodID,index);
            matrixVectorProductRandom(Jacobian3,Rx,Gk_cur3,MethodID,index);
        }

        fk_cur=fk_new;

        BFGS_step(&fk_new,Rx,Xk_cur1,Gk_cur1,Hk_cur1,Jacobian1,Alpha,Mu,Sigma,ImgWrapped,ImgWrappedGrad,Nnum,1,MethodID,index);
        BFGS_step(&fk_new,Rx,Xk_cur2,Gk_cur2,Hk_cur2,Jacobian2,Alpha,Mu,Sigma,ImgWrapped,ImgWrappedGrad,2*Nnum,2,MethodID,index);
        BFGS_step(&fk_new,Rx,Xk_cur3,Gk_cur3,Hk_cur3,Jacobian3,Alpha,Mu,Sigma,ImgWrapped,ImgWrappedGrad,Nnum,3,MethodID,index);

        if (fk_new>fk_cur) break;

        Norm=vectorNorm(Gk_cur1)+vectorNorm(Gk_cur2)+vectorNorm(Gk_cur3);
        printf("Iteracion %d: %lf Norma: %lf\n",k,fk_new,Norm);

        k++;
        string filename="Result"+to_string(k)+".pgm";
        UnWrapImage(Alpha,Mu,Sigma,ImageX);
        NormalizeImage(ImageX,ImageX,255);
        writePGM(&filename[0],ImageX);

        if (fabs(fk_cur-fk_new)<5.0) break;
    }

    FreeMatrixMem(ImageX);
    FreeVectorMem(Alpha);
    FreeVectorMem(Sigma);
    FreeVectorMem(Mu);
    FreeVectorMem(Rx);
    FreeMatrixMem(Jacobian1);
    FreeMatrixMem(Hk_cur1);
    FreeVectorMem(Gk_cur1);
    FreeVectorMem(Xk_cur1);
    FreeMatrixMem(Jacobian2);
    FreeMatrixMem(Hk_cur2);
    FreeVectorMem(Gk_cur2);
    FreeVectorMem(Xk_cur2);
    FreeMatrixMem(Jacobian3);
    FreeMatrixMem(Hk_cur3);
    FreeVectorMem(Gk_cur3);
    FreeVectorMem(Xk_cur3);
    FreeMatrixMem(ImgWrappedGrad[0]);
    FreeMatrixMem(ImgWrappedGrad[1]);
    free(ImgWrappedGrad);
}


#endif // BFGS_H_INCLUDED
