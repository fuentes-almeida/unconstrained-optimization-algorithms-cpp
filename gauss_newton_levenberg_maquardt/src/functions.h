#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

using namespace std;

void CalculateWrappingGradients(Matrix *ImgWrapped,Matrix **ImgWrappedGrad)
{
    int N=ImgWrapped->rows, M=ImgWrapped->cols;

    Matrix *Normalized=RequestMatrixDoubleMem(N,M);
    ImgWrappedGrad[0]=RequestMatrixDoubleMem(N,M);
    ImgWrappedGrad[1]=RequestMatrixDoubleMem(N,M);
    ImageGradient(ImgWrapped,ImgWrappedGrad);

    NormalizeImage(ImgWrappedGrad[0],Normalized,255);
    writePGM("GradX.pgm",Normalized);
    NormalizeImage(ImgWrappedGrad[1],Normalized,255);
    writePGM("GradY.pgm",Normalized);

    for (int i=0;i<N*M;i++)
    {
        double aux1=ImgWrappedGrad[0]->data[i];
        double aux2=ImgWrappedGrad[1]->data[i];

        ImgWrappedGrad[0]->data[i]=atan2(sin(aux1),cos(aux1));
        ImgWrappedGrad[1]->data[i]=atan2(sin(aux2),cos(aux2));
    }

    NormalizeImage(ImgWrappedGrad[0],Normalized,255);
    writePGM("GradXW.pgm",Normalized);
    NormalizeImage(ImgWrappedGrad[1],Normalized,255);
    writePGM("GradYW.pgm",Normalized);

}

double Phi(int i, Vector *Mu, Vector *Sigma, double px, double py)
{
    double dx=px-Mu->data[i], dy=py-Mu->data[Nnum+i];
    double dist=dx*dx+dy*dy;
    double varz=Sigma->data[i]*Sigma->data[i];
    double Phi=exp(-dist/(2*varz));

    return Phi;
}

double VarPhi_1(int i, Vector *Mu, Vector *Sigma, int px, int py)
{
    double dx=px-Mu->data[i], dy=py-Mu->data[Nnum+i];
    double dist=dx*dx+dy*dy;
    double varz=Sigma->data[i]*Sigma->data[i];
    double Phi1=exp(-dist/(2*varz));
    double Phi2=-dx/varz;

    return Phi1*Phi2;
}

double VarPhi_2(int i, Vector *Mu, Vector *Sigma, int px, int py)
{
    double dx=px-Mu->data[i], dy=py-Mu->data[Nnum+i];
    double dist=dx*dx+dy*dy;
    double varz=Sigma->data[i]*Sigma->data[i];
    double Phi1=exp(-dist/(2*varz));
    double Phi2=-dy/varz;

    return Phi1*Phi2;
}

double ObjectiveFunction(Vector *Rx, Matrix **ImgWrappedGrad, Vector *Alpha, Vector *Mu, Vector *Sigma)
{
    int N=imgsize, M=imgsize;
    double result=0;

    for (int px=0;px<N;px++)
        for (int py=0;py<M;py++)
        {
            double sum1=0;
            double sum2=0;
            for (int i=0;i<Nnum;i++)
            {
                sum1+=Alpha->data[i]*VarPhi_1(i,Mu,Sigma,px,py);
                sum2+=Alpha->data[i]*VarPhi_2(i,Mu,Sigma,px,py);
            }
            Rx->data[px*M+py]=ImgWrappedGrad[0]->data[px*M+py]-sum1;
            Rx->data[N*M+px*M+py]=ImgWrappedGrad[1]->data[px*M+py]-sum2;
        }
    result=0.5*vectorInnerProduct(Rx,Rx);

    return result;
}

void ConstructJacobian(Vector **Jacobian, Vector *Alpha, Vector *Mu, Vector *Sigma, int flag)
{
    int N=imgsize, M=imgsize;

    if (flag==1)
    {
        for (int i=0;i<Nnum;i++)
        {
            for (int px=0;px<N;px++)
                for (int py=0;py<M;py++)
                {
                    Jacobian[i]->data[px*M+py]=-VarPhi_1(i,Mu,Sigma,px,py);
                    Jacobian[i]->data[N*M+px*M+py]=-VarPhi_2(i,Mu,Sigma,px,py);
                }
        }
    }
    if (flag==2)
    {
        for (int i=0;i<Nnum;i++)
        {
            for (int px=0;px<N;px++)
                for (int py=0;py<M;py++)
                {
                    double dx=px-Mu->data[i];
                    double dy=py-Mu->data[Nnum+i];
                    double varz=1/(Sigma->data[i]*Sigma->data[i]);

                    Jacobian[i]->data[px*M+py]          =-Alpha->data[i]*varz*(1-varz*dx*dx)*Phi(i,Mu,Sigma,px,py);
                    Jacobian[i]->data[N*M+px*M+py]      = Alpha->data[i]*varz*varz*dx*dy*Phi(i,Mu,Sigma,px,py);
                    Jacobian[Nnum+i]->data[px*M+py]     = Alpha->data[i]*varz*varz*dx*dy*Phi(i,Mu,Sigma,px,py);
                    Jacobian[Nnum+i]->data[N*M+px*M+py] =-Alpha->data[i]*varz*(1-varz*dy*dy)*Phi(i,Mu,Sigma,px,py);
                }
        }
    }
    if (flag==3)
    {
        for (int i=0;i<Nnum;i++)
        {
            for (int px=0;px<N;px++)
                for (int py=0;py<M;py++)
                {
                    double dx=px-Mu->data[i];
                    double dy=py-Mu->data[Nnum+i];
                    double dist=dx*dx+dy*dy;
                    double varz=1/Sigma->data[i];

                    Jacobian[i]->data[px*M+py]=-Alpha->data[i]*Phi(i,Mu,Sigma,px,py)*varz*varz*varz*dx*(2-dist*varz*varz);
                    Jacobian[i]->data[N*M+px*M+py]=-Alpha->data[i]*Phi(i,Mu,Sigma,px,py)*varz*varz*varz*dy*(2-dist*varz*varz);
                }
        }
    }
}

void CalculateGaussianMeans(Vector *Mu, int N, int M)
{
    double delta=imgsize/Num;


    for (int i=0;i<Num;i++)
        for (int j=0;j<Num;j++)
        {
            Mu->data[i*Num+j]=delta/2+i*delta;
            Mu->data[Nnum+i*Num+j]=delta/2+j*delta;
        }
}

void UnWrapImage(Vector *Alpha, Vector *Mu, Vector *Sigma, Matrix *Image)
{
    int N=imgsize, M=imgsize;

    for (int px=0;px<N;px++)
        for (int py=0;py<M;py++)
        {
            double sum=0;
            for (int i=0;i<Nnum;i++)
                sum+=Alpha->data[i]*Phi(i,Mu,Sigma,px,py);

            Image->data[px*M+py]=sum;
        }
}


int SuffDescent(double alphak, double f_xk, Vector *xk, Vector *pk, Vector *gradk, Matrix *ImgWrapped, Matrix **ImgWrappedGrad, Vector *Alpha, Vector *Mu, Vector *Sigma, int size, int flag)
{
    int N=ImgWrapped->rows, M=ImgWrapped->cols;
    double f_xnew, c1=1e-4; //c2=0.9

    double grad_x_p=vectorInnerProduct(gradk,pk);

    Vector *alpha_pk=RequestVectorDoubleMem(size);
    vectorScalarProduct(alpha_pk,pk,alphak);

    Vector *xk_new=RequestVectorDoubleMem(size);
    vectorAdd(xk,alpha_pk,xk_new);

    Vector *R_xnew=RequestVectorDoubleMem(2*N*M);
    if (flag==1)
        f_xnew=ObjectiveFunction(R_xnew,ImgWrappedGrad,xk_new,Mu,Sigma);
    if (flag==2)
        f_xnew=ObjectiveFunction(R_xnew,ImgWrappedGrad,Alpha,xk_new,Sigma);
    if (flag==3)
        f_xnew=ObjectiveFunction(R_xnew,ImgWrappedGrad,Alpha,Mu,xk_new);

    double aux1=f_xk+c1*alphak*grad_x_p;

    FreeVectorMem(alpha_pk);
    FreeVectorMem(xk_new);

    return (f_xnew<=aux1);
}

double LineSearch(double alphak, double f_xk, Vector **Jacobian, Vector *xk, Vector *pk, Matrix *ImgWrapped, Matrix **ImgWrappedGrad, Vector *Alpha, Vector* Mu, Vector *Sigma, int size, int flag)
{
    int N=Jacobian[0]->size;
    Vector *Grad=RequestVectorDoubleMem(size);

    ///se construye el gradiente de la funcion objetivo a partir del Jacobiano
    for (int i=0;i<size;i++)
        for (int j=0;j<N;j++)
            Grad->data[i]+=Jacobian[i]->data[j];

    int maxite=10,count=0,tracking=0,x1=2,x2=5,x3=3;
    double alphaB=alphak;
    while(SuffDescent(alphaB,f_xk,xk,pk,Grad,ImgWrapped,ImgWrappedGrad,Alpha,Mu,Sigma,size,flag)==0 && count<maxite)
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

    if (count==10)
        return alphak;
    else
        return alphaB;
}


void NonLinearLeastSquare(Matrix *ImgWrapped, int MethodID)
{
    int N=imgsize, M=imgsize, k=0;
    double AlphaStepSize=0.8, MuStepSize=0.8, SigmaStepSize=0.8, f_xk;

    Matrix *ImageX=RequestMatrixDoubleMem(N,M);
    Vector *Rx=RequestVectorDoubleMem(2*N*M);
    Vector *Alpha=RequestVectorDoubleMem(Nnum); vectorFill(Alpha,-1);
    Vector *Mu=RequestVectorDoubleMem(2*Nnum);  CalculateGaussianMeans(Mu,N,M);
    Vector *Sigma=RequestVectorDoubleMem(Nnum); vectorFill(Sigma,50);

    Vector **Jacobian=(Vector**)malloc(Nnum*sizeof(Vector*));
    for (int i=0;i<Nnum;i++)
        Jacobian[i]=RequestVectorDoubleMem(2*N*M);

    Vector **JacobianX2=(Vector**)malloc(2*Nnum*sizeof(Vector*));
    for (int i=0;i<2*Nnum;i++)
        JacobianX2[i]=RequestVectorDoubleMem(2*N*M);

    Matrix *JtJ=RequestMatrixDoubleMem(Nnum,Nnum), *JtJX2=RequestMatrixDoubleMem(2*Nnum,2*Nnum);
    Vector *JR=RequestVectorDoubleMem(Nnum),       *JRX2=RequestVectorDoubleMem(2*Nnum);
    Vector *Pk=RequestVectorDoubleMem(Nnum),       *PkX2=RequestVectorDoubleMem(2*Nnum);

    Matrix **ImgWrappedGrad=(Matrix**)malloc(2*sizeof(Matrix*));
    CalculateWrappingGradients(ImgWrapped,ImgWrappedGrad);

    UnWrapImage(Alpha,Mu,Sigma,ImageX);
    NormalizeImage(ImageX,ImageX,255);
    writePGM("Result0.pgm",ImageX);

    while (k<100)
    {
        k++;
            printf("Iteracion %d\n",k);
        ///Optimizacion de Alpha
            f_xk=ObjectiveFunction(Rx,ImgWrappedGrad,Alpha,Mu,Sigma);
            printf("Fitness: %lf\n",f_xk);

            ConstructJacobian(Jacobian,Alpha,Mu,Sigma,1);

            if (MethodID==1)
                GaussNewtonMethod(Jacobian,JtJ,Rx,JR,Pk,Nnum);
            if (MethodID==2)
                LevenbergMaquardtMethod(Jacobian,JtJ,Rx,JR,Pk,Nnum);

            printf("alfa Pk: ");
            for (int i=0;i<Nnum;i++)
                printf("%lf ",Pk->data[i]);
            printf("Norma: %lf\n",vectorNorm(Pk));

            AlphaStepSize=LineSearch(AlphaStepSize,f_xk,Jacobian,Alpha,Pk,ImgWrapped,ImgWrappedGrad,Alpha,Mu,Sigma,Nnum,1);

            vectorScalarProduct(Pk,Pk,AlphaStepSize);
            vectorAdd(Pk,Alpha,Alpha);

        ///Optimizacion de Mu
            f_xk=ObjectiveFunction(Rx,ImgWrappedGrad,Alpha,Mu,Sigma);
            printf("Fitness: %lf\n",f_xk);

            ConstructJacobian(JacobianX2,Alpha,Mu,Sigma,2);

            if (MethodID==1)
                GaussNewtonMethod(JacobianX2,JtJX2,Rx,JRX2,PkX2,2*Nnum);
            if (MethodID==2)
                LevenbergMaquardtMethod(JacobianX2,JtJX2,Rx,JRX2,PkX2,2*Nnum);

            printf("Mu Pk: ");
            for (int i=0;i<Nnum;i++)
                printf("%lf,%lf ",PkX2->data[i],PkX2->data[Nnum+i]);
            printf("Norma: %lf\n",vectorNorm(PkX2));

            MuStepSize=LineSearch(MuStepSize,f_xk,JacobianX2,Mu,PkX2,ImgWrapped,ImgWrappedGrad,Alpha,Mu,Sigma,2*Nnum,2);

            //vectorScalarProduct(PkX2,PkX2,MuStepSize);
            vectorAdd(PkX2,Mu,Mu);

        ///Optimizacion de Sigma
            f_xk=ObjectiveFunction(Rx,ImgWrappedGrad,Alpha,Mu,Sigma);
            printf("Fitness: %lf\n",f_xk);

            ConstructJacobian(Jacobian,Alpha,Mu,Sigma,3);

            if (MethodID==1)
                GaussNewtonMethod(Jacobian,JtJ,Rx,JR,Pk,Nnum);
            if (MethodID==2)
                LevenbergMaquardtMethod(Jacobian,JtJ,Rx,JR,Pk,Nnum);

            printf("Sigma Pk: ");
            for (int i=0;i<Nnum;i++)
                printf("%lf ",Pk->data[i]);
            printf("Norma: %lf\n\n",vectorNorm(Pk));

            SigmaStepSize=LineSearch(SigmaStepSize,f_xk,Jacobian,Sigma,Pk,ImgWrapped,ImgWrappedGrad,Alpha,Mu,Sigma,Nnum,3);

            //vectorScalarProduct(Pk,Pk,SigmaStepSize);
            vectorAdd(Pk,Sigma,Sigma);

        string filename="Result"+to_string(k)+".pgm";

        UnWrapImage(Alpha,Mu,Sigma,ImageX);
        NormalizeImage(ImageX,ImageX,255);
        writePGM(&filename[0],ImageX);
    }

    FreeMatrixMem(ImageX);
    FreeVectorMem(Alpha);
    FreeVectorMem(Sigma);
    FreeVectorMem(Mu);
    FreeVectorMem(Rx);

    for (int i=0;i<Nnum;i++)
        FreeVectorMem(Jacobian[i]);
    for (int i=0;i<2*Nnum;i++)
        FreeVectorMem(JacobianX2[i]);

    free(Jacobian);
    free(JacobianX2);

    FreeMatrixMem(JtJ);
    FreeMatrixMem(JtJX2);
    FreeVectorMem(JR);
    FreeVectorMem(JRX2);
    FreeVectorMem(Pk);
    FreeVectorMem(PkX2);

    FreeMatrixMem(ImgWrappedGrad[0]);
    FreeMatrixMem(ImgWrappedGrad[1]);
    free(ImgWrappedGrad);
}

#endif // FUNCTIONS_H_INCLUDED
