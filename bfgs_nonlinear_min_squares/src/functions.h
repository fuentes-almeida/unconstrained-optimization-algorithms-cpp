#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

using namespace std;

void RandomVector(int rnum, int rmax, int *index)
{
    for (int i=0;i<rnum;i++)
        index[i]=0;

    for (int i=0;i<rnum;i++)
    {
        if (1.0*rand()/RAND_MAX<0.5)
            index[i]=2*i;
        else
            index[i]=2*i+1;
    }
}

void matrixVectorProductRandom(Matrix * A, Vector * x, Vector * b, int methodID, int *index)
{
    if (methodID==1)
    {
        int N=A->rows, M=A->cols;
        for (int i=0;i<N;i++)
        {
            double sum=0;
            for(int j=0;j<M;j++)
                sum+=A->data[i*M+j]*x->data[j];
            b->data[i]=sum;
        }
    }
    if (methodID==2)
    {
        int N=A->rows, M=A->cols;
        for (int i=0;i<N;i++)
        {
            double sum=0;

            for(int j=0;j<M;j++)
                sum+=A->data[i*M+j]*x->data[index[j%(M/2)]];

            b->data[i]=sum;
        }
    }
}

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

void ConstructJacobian(Matrix *Jacobian, Vector *Alpha, Vector *Mu, Vector *Sigma, int flag, int methodID, int *index)
{
    int N=imgsize, M=imgsize;
    int dim=M*N, cols=2*dim;

    if (methodID==1)
    {
        if (flag==1)
        {
            for (int i=0;i<Nnum;i++)
            {
                for (int px=0;px<N;px++)
                    for (int py=0;py<M;py++)
                    {
                        Jacobian->data[i*cols+px*M+py]     =-VarPhi_1(i,Mu,Sigma,px,py);
                        Jacobian->data[i*cols+dim+px*M+py] =-VarPhi_2(i,Mu,Sigma,px,py);
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

                        Jacobian->data[i*cols+px*M+py]            =-Alpha->data[i]*varz*(1-varz*dx*dx)*Phi(i,Mu,Sigma,px,py);
                        Jacobian->data[i*cols+dim+px*M+py]        = Alpha->data[i]*varz*varz*dx*dy*Phi(i,Mu,Sigma,px,py);
                        Jacobian->data[(i+Nnum)*cols+px*M+py]     = Alpha->data[i]*varz*varz*dx*dy*Phi(i,Mu,Sigma,px,py);
                        Jacobian->data[(i+Nnum)*cols+dim+px*M+py] =-Alpha->data[i]*varz*(1-varz*dy*dy)*Phi(i,Mu,Sigma,px,py);
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

                        Jacobian->data[i*cols+px*M+py]     =-Alpha->data[i]*Phi(i,Mu,Sigma,px,py)*varz*varz*varz*dx*(2-dist*varz*varz);
                        Jacobian->data[i*cols+dim+px*M+py] =-Alpha->data[i]*Phi(i,Mu,Sigma,px,py)*varz*varz*varz*dy*(2-dist*varz*varz);
                    }
            }
        }
    }

    if (methodID==2)
    {
        if (flag==1)
        {
            for (int i=0;i<Nnum;i++)
            {
                for (int j=0;j<dim/2;j++)
                    {
                        int px=index[j]/N;
                        int py=index[j]%M;
                        Jacobian->data[i*dim+j]       =-VarPhi_1(i,Mu,Sigma,px,py);
                        Jacobian->data[i*dim+dim/2+j] =-VarPhi_2(i,Mu,Sigma,px,py);
                    }
            }
        }
        if (flag==2)
        {
            for (int i=0;i<Nnum;i++)
            {
                for (int j=0;j<dim/2;j++)
                    {
                        int px=index[j]/N;
                        int py=index[j]%M;

                        double dx=px-Mu->data[i];
                        double dy=py-Mu->data[Nnum+i];
                        double varz=1/(Sigma->data[i]*Sigma->data[i]);

                        Jacobian->data[i*dim+j]              =-Alpha->data[i]*varz*(1-varz*dx*dx)*Phi(i,Mu,Sigma,px,py);
                        Jacobian->data[i*dim+dim/2+j]        = Alpha->data[i]*varz*varz*dx*dy*Phi(i,Mu,Sigma,px,py);
                        Jacobian->data[(i+Nnum)*dim+j]       = Alpha->data[i]*varz*varz*dx*dy*Phi(i,Mu,Sigma,px,py);
                        Jacobian->data[(i+Nnum)*dim+dim/2+j] =-Alpha->data[i]*varz*(1-varz*dy*dy)*Phi(i,Mu,Sigma,px,py);
                    }
            }
        }
        if (flag==3)
        {
            for (int i=0;i<Nnum;i++)
            {
                for (int j=0;j<dim/2;j++)
                    {
                        int px=index[j]/N;
                        int py=index[j]%M;

                        double dx=px-Mu->data[i];
                        double dy=py-Mu->data[Nnum+i];
                        double dist=dx*dx+dy*dy;
                        double varz=1/Sigma->data[i];

                        Jacobian->data[i*dim+j]       =-Alpha->data[i]*Phi(i,Mu,Sigma,px,py)*varz*varz*varz*dx*(2-dist*varz*varz);
                        Jacobian->data[i*dim+dim/2+j] =-Alpha->data[i]*Phi(i,Mu,Sigma,px,py)*varz*varz*varz*dy*(2-dist*varz*varz);
                    }
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

int SmoothWolfe(double alphak, double fk, Vector *Xk, Vector *Gk, Vector *Pk, Matrix *ImgWrapped, Matrix **ImgWrappedGrad, Vector *Alpha, Vector *Mu, Vector *Sigma, int size, int flag, int methodID, int *index)
{
    int N=ImgWrapped->rows, M=ImgWrapped->cols, tam;
    double fk_new, c1=1e-4, c2=0.9;

    if (methodID==1)
        tam=2*N*M;
    if (methodID==2)
        tam=N*M;

    Vector *alpha_p=RequestVectorDoubleMem(size);
    Vector *x_new=RequestVectorDoubleMem(size);
    Vector *grad_xnew=RequestVectorDoubleMem(size);

    vectorScalarProduct(alpha_p,Pk,alphak);
    vectorAdd(Xk,alpha_p,x_new);

    double grad_x_p=vectorInnerProduct(Gk,Pk);

    Vector *R_xnew=RequestVectorDoubleMem(2*N*M);
    if (flag==1)
        fk_new=ObjectiveFunction(R_xnew,ImgWrappedGrad,x_new,Mu,Sigma);
    if (flag==2)
        fk_new=ObjectiveFunction(R_xnew,ImgWrappedGrad,Alpha,x_new,Sigma);
    if (flag==3)
        fk_new=ObjectiveFunction(R_xnew,ImgWrappedGrad,Alpha,Mu,x_new);

    double aux1=fk+c1*alphak*grad_x_p;

    ///nuevo gradiente
    Matrix *Jacobian=RequestMatrixDoubleMem(size,tam);
    if (flag==1)
        ConstructJacobian(Jacobian,x_new,Mu,Sigma,flag,methodID,index);
    if (flag==2)
        ConstructJacobian(Jacobian,Alpha,x_new,Sigma,flag,methodID,index);
    if (flag==3)
        ConstructJacobian(Jacobian,Alpha,Mu,x_new,flag,methodID,index);

    matrixVectorProductRandom(Jacobian,R_xnew,grad_xnew,methodID,index);
    double aux2=vectorInnerProduct(grad_xnew,Pk);

    FreeVectorMem(alpha_p);
    FreeVectorMem(x_new);
    FreeVectorMem(grad_xnew);
    FreeVectorMem(R_xnew);
    FreeMatrixMem(Jacobian);

    return (fk_new<=aux1 && aux2>=c2*grad_x_p);
}

double LineSearch(double alphak, double fk, Vector *Xk, Vector *Gk, Vector *Pk, Matrix *ImgWrapped, Matrix **ImgWrappedGrad, Vector *Alpha, Vector* Mu, Vector *Sigma, int size, int flag, int methodID, int *index)
{
    int maxite=10,count=0,tracking=0,x1=2,x2=5,x3=3;
    double alphaB=alphak;
    while(SmoothWolfe(alphaB,fk,Xk,Gk,Pk,ImgWrapped,ImgWrappedGrad,Alpha,Mu,Sigma,size,flag,methodID,index)==0 && count<maxite)
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


#endif // FUNCTIONS_H_INCLUDED
