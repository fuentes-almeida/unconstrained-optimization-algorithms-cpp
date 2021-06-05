#ifndef FUNCTION_H_INCLUDED
#define FUNCTION_H_INCLUDED

///funcion que calcula el peso sobre los vecinos de un pixel
double weight(double Urs, double k, int wflag)
{
    double ws=1.0;
    if (wflag==0)
        ws=k/(k+Urs*Urs);

    double result=ws*ws*Urs;
    return result;
}

///funcion que calcula el potencial de acuerdo a la expresion dada en la especificacion
double potential(double Urs, double k, int wflag)
{
    double ws=1.0;
    if (wflag==0)
        ws=k/(k+Urs*Urs);
    double result=ws*ws*Urs*Urs+(1-ws)*(1-ws)*k;
    return result;
}

///Funcion que calcula el valor de la funcion del modelo a partir de una imagen obtenida
double ModelFunction(Matrix *ImageX, Matrix *ImageY, double lambda, double k, int wflag)
{
    int N=ImageX->rows;
    int M=ImageX->cols;

    double error=0.0;

    for(int i=0; i<N; i++)
        for (int j=0; j<M; j++)
        {
            error+=0.5*pow(ImageX->data[i*M+j]-ImageY->data[i*M+j],2);
            ///se utiliza la definicion de vecindad especificada en la tarea
            if (i>0)    error+=0.5*lambda*potential(ImageX->data[i*M+j]-ImageX->data[(i-1)*M+j],k,wflag);
            if (i<N-1)  error+=0.5*lambda*potential(ImageX->data[i*M+j]-ImageX->data[(i+1)*M+j],k,wflag);
            if (j>0)    error+=0.5*lambda*potential(ImageX->data[i*M+j]-ImageX->data[i*M+(j-1)],k,wflag);
            if (j<M-1)  error+=0.5*lambda*potential(ImageX->data[i*M+j]-ImageX->data[i*M+(j+1)],k,wflag);
        }
    return error;
}

///Gradiente del modelo con respecto a cada pixel de la imagen X
void ModelFunctionGradient(Matrix *ImageX, Matrix *ImageY, Matrix *Grad, double lambda, double k, int wflag)
{
    int N=ImageX->rows;
    int M=ImageX->cols;

    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
        {
            double aux=ImageX->data[i*M+j]-ImageY->data[i*M+j];
            ///se utiliza la definicion de vecindad especificada en la tarea
            if (i>0)    aux+=lambda*weight(ImageX->data[i*M+j]-ImageX->data[(i-1)*M+j],k,wflag);
            if (i<N-1)  aux+=lambda*weight(ImageX->data[i*M+j]-ImageX->data[(i+1)*M+j],k,wflag);
            if (j>0)    aux+=lambda*weight(ImageX->data[i*M+j]-ImageX->data[i*M+(j-1)],k,wflag);
            if (j<M-1)  aux+=lambda*weight(ImageX->data[i*M+j]-ImageX->data[i*M+(j+1)],k,wflag);
            Grad->data[i*M+j]=aux;
        }
}

#endif // FUNCTION_H_INCLUDED
