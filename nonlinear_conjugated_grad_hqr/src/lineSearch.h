#ifndef LINESEARCH_H_INCLUDED
#define LINESEARCH_H_INCLUDED

///condicion fuerte de Wolfe para el calculo del tamanio de paso alpha
int StrongWolfe(double alpha, double gkpk, Matrix *x, Matrix *y, Matrix *p, Matrix *grad, double fx, double lambda, double k)
{
    double c1=0.01,c2=0.2;
    int N=x->rows, M=x->cols;
    Matrix *alpha_p=RequestMatrixDoubleMem(N,M);
    Matrix *x_new=RequestMatrixDoubleMem(N,M);
    Matrix *grad_xnew=RequestMatrixDoubleMem(N,M);

    matrixScalarProduct(alpha_p,p,alpha);
    matrixAdd(x,alpha_p,x_new);

    double aux1=fx+c1*alpha*gkpk;

    ModelFunctionGradient(x_new,y,grad_xnew,lambda,k,0);

    double aux2=matrixInnerProduct(grad_xnew,p);

    FreeMatrixMem(alpha_p);
    FreeMatrixMem(x_new);
    FreeMatrixMem(grad_xnew);

    double fx_new=ModelFunction(x_new,y,lambda,k,0);

    return (fx_new<=aux1 && fabs(aux2)<=c2*fabs(gkpk));
}

///modelo cuadratico para el calculo del tamanio de paso alpha
double CuadraticModel(double alpha1, Matrix *x, Matrix *y, Matrix *p, Matrix *grad, double fx, double lambda, double k)
{
    int N=x->rows, M=x->cols;
    double gkpk=matrixInnerProduct(grad,p);

    ///si el valor de alpha recibido como parametro ya cumple la condicion se devuelve
    if (StrongWolfe(alpha1,gkpk,x,y,p,grad,fx,lambda,k)==1)
        return alpha1;

    ///si no cumple, se calcula una alpha a partir de un modelo cuadratico
        double phi01=fx;
        double phi02=matrixInnerProduct(grad,p);

        Matrix *alpha_p=RequestMatrixDoubleMem(N,M);
        Matrix *x_new=RequestMatrixDoubleMem(N,M);

        matrixScalarProduct(alpha_p,p,alpha1);
        matrixAdd(x,alpha_p,x_new);

        double phi_alpha1=ModelFunction(x_new,y,lambda,k,0);
        double alpha2=(alpha1*alpha1*phi02)/(2.0*(phi01+phi02*alpha1-phi_alpha1));

    ///si el valor de alpha calculado ya cumple la condicion se devuelve
    if (StrongWolfe(alpha2,gkpk,x,y,p,grad,fx,lambda,k)==1)
        return alpha2;

    FreeMatrixMem(alpha_p);
    FreeMatrixMem(x_new);

    ///si ningun valor cumplio la condicion se regresa el valor de alpha recibido
    return alpha1;
}

#endif // LINESEARCH_H_INCLUDED
