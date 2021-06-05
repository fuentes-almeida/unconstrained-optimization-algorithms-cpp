#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "MatrixInterface.h"
#include "CommonTypes.h"

Vector *RequestVectorDoubleMem(int data_size)
{
    Vector *V=(Vector*)malloc(sizeof(Vector));
    V->size=data_size;
    V->data=(double*)calloc(data_size,sizeof(double));
    return V;
}

Matrix *RequestMatrixDoubleMem(int data_size1, int data_size2)
{
    Matrix *M=(Matrix*)malloc(sizeof(Matrix));
    M->rows=data_size1;
    M->cols=data_size2;
    M->size=data_size1*data_size2;
    M->data=(double*)calloc(data_size1*data_size2,sizeof(double));
    return M;
}

void FreeVectorMem(Vector *V)
{
    free(V->data);
    free(V);
}

void FreeMatrixMem(Matrix *M)
{
    free(M->data);
    free(M);
}

/***************************************************************************Operaciones entre vectores y matrices******************/

/*!
\brief Se encarga de realizar la operaci\'on vectorial c = a - b
\param a Vector de entrada
\param b Vector de entrada
\param c Vector de resultado
*/
void vectorSubstract(Vector *a, Vector *b, Vector *c)
{
    int dim=c->size;
    for (int i=0;i<dim;i++)
        c->data[i]=a->data[i]-b->data[i];
}

/*!
\brief Se encarga de realizar la operaci\'on vectorial c = a + b
\param a Vector de entrada
\param b Vector de entrada
\param c Vector de resultado
*/
void vectorAdd(Vector *a, Vector *b, Vector *c)
{
    int dim=c->size;
    for (int i=0;i<dim;i++)
        c->data[i]=a->data[i]+b->data[i];
}


/*!
\brief Se encarga de realizar la operaci\'on vectorial a = b*s
\param a Vector de entrada
\param b Vector de resultado
\param s Scalar
*/
void vectorScalarProduct(Vector *a, Vector * b, double s)
{
    int dim=a->size;
    for (int i=0;i<dim;i++)
        a->data[i]=s*b->data[i];
}


/*!
\brief Se encarga de calcular la distancia entre dos vector dist  = || a - b|| ^{1/2}
\param a Vector de entrada
\param b Vector de resultado
\return d Disctancia
*/
double vectorDistance(Vector *a, Vector *b)
{
	return 0;
}


/*!
\brief Se encarga de calcular la norma al cuadrado del vector, norm = ||a||^2
\param a Vector de entrada.
\return n Norma del vector.
*/
double vectorNorm2(Vector *a)
{
    double norm=0;
    for (int i=0;i<a->size;i++)
        norm+=a->data[i]*a->data[i];
	return norm;
}


/*!
\brief Se encarga de calcular la norma  del vector, norm = ||a||
\param a Vector de entrada.
\return n Norma del vector.
*/
double vectorNorm(Vector *a)
{
    double norm=0;
    for (int i=0;i<a->size;i++)
        norm+=a->data[i]*a->data[i];
	return sqrt(norm);
}

///Calcula la norma de una matriz como si fuera un vector
double matrixNorm(Matrix *M)
{
    double norm=0;
    int dim=M->size;
    for (int i=0;i<dim;i++)
        norm+=M->data[i]*M->data[i];
    return sqrt(norm);
}

/*!
\brief Se encarga de calcular el vector unitario de a, b = a / ||a||.
\param a Vector a calcular su vector unitario.
\param b Vector unitario resultante.
*/
void vectorNormalize(Vector *a, Vector * b)
{
}


/*!
\brief se encarga de calcular el producto punto entre dos vectores, c = a^T * b.
\param a Primer vector.
\param b Segundo vector.
\return c Valor resultante del producto punto.
*/
double vectorInnerProduct(Vector * a, Vector * b)
{
    int dim=a->size;
    double result=0;
    for (int i=0;i<dim;i++)
        result+=a->data[i]*b->data[i];

	return result;
}

///Calcula el producto interno de dos matrices como si fueran vectores
double matrixInnerProduct(Matrix *A, Matrix *B)
{
    int dim=A->size;
    double result=0.0;
    for (int i=0;i<dim;i++)
        result+=A->data[i]*B->data[i];
	return result;
}


/*!
\brief se encarga de calcular el producto externo entre dos vectores, C = a * b^T.
\param a Primer vector.
\param b Segundo vector.
\param C Matriz resultante el producto externo.
*/
void vectorOuterProduct(Vector * a, Vector * b, Matrix *C)
{
    for (int i=0;i<C->rows;i++)
        for (int j=0;j<C->cols;j++)
        {
            C->data[i*C->cols+j]=a->data[i]*b->data[j];
        }
}


/*!
\brief se encarga de encontrar el m\'inimo elemento del vector a.
\param a Vector a calcular el m\'inimo.
\return b M\'inimo elemento del vector.
*/
double vectorMin(Vector *a)
{
	return 0;
}


/*!
\brief se encarga de encontrar el m\'aximo elemento del vector a.
\param a Vector a calcular el m\'aximo.
\return b M\'aximo elemento del vector.
*/
double vectorMax(Vector *a)
{
	return 0;
}


/*!
\brief Se encarga de realizar la operaci\'on matricial C = A + B
\param A Matriz de entrada
\param B Matriz de entrada
\param C Matriz resultante
*/
void matrixAdd(Matrix * A, Matrix * B, Matrix * C)
{
    for (int i=0;i<C->rows;i++)
        for (int j=0;j<C->cols;j++)
            C->data[i*C->cols+j]=A->data[i*A->cols+j]+B->data[i*B->cols+j];
}


/*!
\brief Se encarga de realizar la operaci\'on matricial C = A - B
\param A Matriz de entrada
\param B Matriz de entrada
\param C Matriz resultante
*/
void matrixSubstract(Matrix * A, Matrix * B, Matrix * C)
{
    for (int i=0;i<C->rows;i++)
        for (int j=0;j<C->cols;j++)
            C->data[i*C->cols+j]=A->data[i*A->cols+j]-B->data[i*B->cols+j];
}


/*!
\brief Se encarga de realizar la operaci\'on matricial A = B*s
\param A Matriz de entrada
\param B Matriz resultante
\param s Scalar
*/
void matrixScalarProduct(Matrix * A, Matrix * B, double s)
{
    for (int i=0;i<A->rows;i++)
        for (int j=0;j<A->cols;j++)
            A->data[i*A->cols+j]=B->data[i*B->cols+j]*s;

}


/*!
\brief Se encarga de realizar la operaci\'on sobre la diagonal de A,  A(i,i) = B(i,i) + s
\param A Matriz de entrada
\param B Matriz resultante
\param s Scalar
*/
void matrixAddToDiagonal(Matrix * A, Matrix * B, double s)
{
    for (int i=0;i<A->rows;i++)
        A->data[i*A->cols+i]=B->data[i*B->cols+i]+s;

}

/*!
\brief Se encarga de realizar la operaci\'on matricial B = A^T
\param a Vector de entrada
\param b Vector de resultado
\param s Scalar
*/
void matrixTranspose(Matrix * A, Matrix * B)
{
    for (int i=0;i<B->rows;i++)
        for (int j=0;j<B->cols;j++)
            B->data[i*B->cols+j]=A->data[j*A->cols+i];
}

/*!
\brief Se encarga de realizar la operaci\'on matricial C = A*B
\param A Matriz de entrada
\param B Matriz de entrada
\param C Matriz resultante
*/
void matrixProduct(Matrix * A, Matrix * B, Matrix * C)
{
    matrixFill(C,0);
    for (int i=0;i<A->rows;i++)
        for (int j=0;j<B->cols;j++)
            for (int k=0;k<B->rows;k++)
                C->data[i*C->cols+j]+=A->data[i*A->cols+k]*B->data[k*B->cols+j];
}


/*!
\brief se encarga de encontrar el m\'inimo elemento de la matriz A.
\param A Matriz a calcular el m\'inimo.
\return b M\'inimo elemento de la matriz.
*/
double matrixMin(Matrix *A)
{
    double mini=1e9;
    for (int i=0;i<A->size;i++)
        if (A->data[i]<mini)
            mini=A->data[i];

	return mini;
}


/*!
\brief se encarga de encontrar el m\'aximo elemento de la matriz A.
\param A Matriz a calcular el m\'aximo.
\return b M\'aximo elemento de la matriz.
*/
double matrixMax(Matrix *A)
{
    double maxi=0;
    for (int i=0;i<A->size;i++)
        if (A->data[i]>maxi)
            maxi=A->data[i];

	return maxi;
}

/*!
\brief se encarga de calcular el producto entre una matriz y un vector, b = A * x.
\param A Matriz de entrada.
\param x Vector de incognitas.
\param b Vector resultante de la multiplicaci\'on.
*/

void matrixVectorProduct(Matrix * A, Vector * x, Vector * b)
{
    int dim=x->size;
    for (int i=0;i<dim;i++){
            double suma=0;
        for(int j=0;j<dim;j++){
            suma+=A->data[i*dim+j]*x->data[j];
        }
        b->data[i]=suma;
    }
}


/*!
\brief se encarga de calcular la forma cuadratica entre una matriz y un vector, c = x^T * A * x.
\param A Matriz de entrada.
\param x Vector de incognitas.
\return c Valor resultante de la forma cuadratica
*/
double quadraticFormProduct(Matrix * A, Vector * x)
{
    Vector *Aux=RequestVectorDoubleMem(x->size);
    matrixVectorProduct(A,x,Aux);
            for (int i=0;i<9;i++)
            printf("%lf ",Aux->data[i]);
        printf("\n");
    double product=vectorInnerProduct(Aux,x);
    FreeVectorMem(Aux);
	return product;
}

/************************************************************Calculo Vectorial******************************************************************************/

double PartialDeriv(Vector *x,int i,double h,funcion fx)
{
    if (i==0)
        return -400*x->data[0]*(x->data[1]-x->data[0]*x->data[0])-2*(1-x->data[0]);
    else if (i==1)
        return 200*(x->data[1]-x->data[0]*x->data[0]);
    else
        return 0;
}

double MixPartialDeriv(Vector *x,int i,int j,double h1, double h2,funcion fx)
{
    if (i==0 && j==0)
        return -400*(x->data[1]-x->data[0]*x->data[0])+800*x->data[0]*x->data[0]+2;
    else if (i==1 && j==1)
        return 200;
    else if ((i==0 && j==1) || (i==1 && j==0) )
        return -400*x->data[0];
    else
        return 0;

}




/***********************************************************Factorizaciones y solvers para sistemas de ecuaciones Ax = b************************************/

//-----------------------METODO DE CROUTH---------------------------------//

void SolveTriSupMatrix(Matrix *MatrixSup,Vector *x,Vector *b,int n)
{
    int i,j;
    double sum=0;
    x->data[n-1]=b->data[n-1];
    for(j=n-2;j>=0;j--){
                sum=0;
        for(i=j+1;i<n;i++){
            sum=sum + MatrixSup->data[j*n+i]*x->data[i];
        }
        x->data[j]=b->data[j]-sum;
    }
}

void SolveTriInfMatrix(Matrix *MatrixInf,Vector *x,Vector *b,int n)
{
    int i,j;
    double sum=0;
    x->data[0]=b->data[0]/MatrixInf->data[0*n+0];
    for(j=1;j<n;j++){
             sum=0;
        for(i=0;i<=j-1;i++){
            sum=sum + MatrixInf->data[j*n+i]*x->data[i];}
        x->data[j]=(b->data[j]-sum)/MatrixInf->data[j*n+j];
    }
}


void SolveCrouthSystem(Matrix *A, Vector *b, Vector *x,int n)
{
    croutDecomposition(A,n);
    SolveTriInfMatrix(A,x,b,n);
    SolveTriSupMatrix(A,x,x,n);
}

/*!\brief Obtiene la factorizaci\'on A = LU de la matriz de entrada A, usando el m\'etodo de Crout.
\param A Matriz a descomponer.
\param L Matriz a inferior.
\param U Matriz a superior.
*/
void croutDecomposition(Matrix * A,int n)
{
    int i,j,k;
    double sum;
    for (i=0;i<n;i++){
    sum=0;
    for(k=0;k<i;k++)
        sum=sum+A->data[i*n+k]*A->data[k*n+i];
    A->data[i*n+i]=A->data[i*n+i]-sum;

    for(j=i+1;j<n;j++){
            sum=0;
        for(k=0;k<i;k++)
            sum=sum+A->data[j*n+k]*A->data[k*n+i];
    A->data[j*n+i]=A->data[j*n+i]-sum;}

    for(j=i+1;j<n;j++){
            sum=0;
        for (k=0;k<i;k++)
            sum=sum+A->data[i*n+k]*A->data[k*n+j];
    A->data[i*n+j]=(A->data[i*n+j]-sum)/A->data[i*n+i];}
    }
}


/*!\brief Obtiene la factorizaci\'on LU de la matriz de entrada, usando el m\'etodo de Dolitle.
\param A Matriz a descomponer.
\param L Matriz a inferior.
\param U Matriz a superior.
*/
void dolitleDecomposition(Matrix * A, Matrix * L, Matrix * U)
{
}


/*!
\brief Obtiene la factorizacion A = LL^T de la matriz bandada A de entrada, usando el metodo de cholesky.
\param A Matriz a descomponer.
\param L Matriz a inferior.
\param bandwidth Dimension de la banda.
*/
void choleskyDecomposition(Matrix * A, Matrix *L, int bandwidth)
{
}


/*!\brief Obtiene la factorizaci\'on LDL^T de la matriz sim\'etrica de entrada, usando el m\'etodo de cholesky.
\param A Matriz a descomponer.
\param L Matriz a inferior.
*/
void choleskyDecomposition(Matrix * A)
{
    int i,j,k;
    double sum=0;
    int n=A->rows;
    for (i=0;i<n;i++)
    {
        sum=0;
        for(k=0;k<i;k++)
            sum=sum+A->data[i*n+k]*A->data[i*n+k]*A->data[k*n+k];
        A->data[i*n+i]=A->data[i*n+i]-sum;

        for(j=i+1;j<n;j++)
        {
            sum=0;
            for(k=0;k<i;k++)
                sum=sum+A->data[i*n+k]*A->data[j*n+k]*A->data[k*n+k];

            A->data[j*n+i]=(A->data[j*n+i]-sum)/A->data[i*n+i];
            A->data[i*n+j]=A->data[j*n+i];
        }
    }
}

/*!
\brief Resuelve un sistema lineal Ax=b, ssi A es una matriz triangular inferior tal que A(i,i) != 0.
\param A Matriz a resolver.
\param b Vector de soluciones.
\param x Vector de incognitas.
*/
void forwardSubstitution(Matrix * A, Vector *b, Vector * x)
{
    int i,j;
    double sum=0;
    int n=x->size;

    x->data[0]=b->data[0];
    for(j=1;j<n;j++){
             sum=0;
        for(i=0;i<=j-1;i++){
            sum=sum + A->data[j*n+i]*x->data[i];}
        x->data[j]=b->data[j]-sum;
    }
}

/*!
\brief Resuelve un sistema lineal Ax=b, ssi A es una matriz triangular superior tal que A(i,i) != 0.
\param A Matriz a resolver.
\param b Vector de soluciones.
\param x Vector de incognitas.
*/
void backwardSubstitution(Matrix * A, Vector *b, Vector * x)
{
    int i,j;
    double sum=0;
    int n=x->size;

    x->data[n-1]=b->data[n-1];
    for(j=n-2;j>=0;j--){
                sum=0;
        for(i=j+1;i<n;i++){
            sum=sum + A->data[j*n+i]*x->data[i];
        }
        x->data[j]=b->data[j]-sum;
    }
}

void SolveDiagMatrix(Matrix *D,Vector *x,Vector *b)
{
    int i;
    int n=x->size;
    for(i=0;i<n;i++){
        x->data[i]=b->data[i]/D->data[i*n+i];}

}


/*!
\brief Resuelve un sistema lineal Ax=b, por el metodo de cholesky ssi A es una Matriz bandada.
\param A Matriz a resolver.
\param b Vector de soluciones.
\param x Vector de incognitas.
\param bandwidth Dimension de la banda.
*/
void choleskySolve(Matrix * A, Vector * b, Vector * x, int bandwidth)
{
}


/*!
\brief Resuelve un sistema lineal Ax=b, por el metodo de cholesky ssi A = LL^T.
\param A Matriz a resolver.
\param b Vector de soluciones.
\param x Vector de incognitas.
*/
void choleskySolve(Matrix * A, Vector * b, Vector * x)
{
    choleskyDecomposition(A);
    forwardSubstitution(A,b,x);
    SolveDiagMatrix(A,x,x);
    backwardSubstitution(A,x,x);
}


/*!\brief Resuelve un sistema lineal Ax=b por medio del m\'etodo de gauss-jordan.
\param A Matriz a resolver.
\param b Vector de soluciones.
\param x Vector de incognitas.
*/
void gaussJordan(Matrix * A, Vector * b, Vector * x)
{
}

///Metodo de Gradiente Conjugado
void ConjugatedGradient(Matrix *A, Vector *b, Vector *x)
{
    Vector *gradk=RequestVectorDoubleMem(x->size);
    Vector *gradk_New=RequestVectorDoubleMem(x->size);
    Vector *pk=RequestVectorDoubleMem(x->size);
    Vector *aux=RequestVectorDoubleMem(x->size);

    vectorSubstract(gradk,b,gradk);
    vectorScalarProduct(pk,gradk,-1);

    double epsilon=1e-9,alpha=0,beta=0;
    int Itmax=10*x->size,k=0;

    while ((vectorNorm(gradk)>epsilon) & (k<Itmax))
    {
        k++;
        double gkTpk=vectorInnerProduct(gradk,pk);
        double pkTApk=quadraticFormProduct(A,pk);
        alpha=-gkTpk/pkTApk;

        vectorScalarProduct(aux,pk,alpha);
        vectorAdd(x,aux,x);

        vectorFill(aux,0);
        matrixVectorProduct(A,pk,aux);
        vectorScalarProduct(aux,aux,alpha);

        vectorAdd(gradk,aux,gradk_New);

        double gkTgk_New=vectorInnerProduct(gradk_New,gradk_New);
        double gkTgk=vectorInnerProduct(gradk,gradk);
        beta=gkTgk_New/gkTgk;

        vectorScalarProduct(pk,pk,beta);
        vectorSubstract(pk,gradk_New,pk);

        vectorCopy(gradk_New,gradk);
    }

}

/***********************************************************************************************Funciones de utileria********************/

double Rosenbrock(Vector *x)
{
    int dim=x->size;
    double fx=0;
    for (int i=0;i<dim-1;i++)
        fx+=100*(x->data[i+1]-x->data[i]*x->data[i])*(x->data[i+1]-x->data[i]*x->data[i])+(1-x->data[i])*(1-x->data[i]);
    return fx;
}

double Euclidean(Vector *v1,Vector *v2)
{
    int dim=v1->size;
    double dist=0;
    for (int i=0;i<dim;i++)
        dist+=(v1->data[i]-v2->data[i])*(v1->data[i]-v2->data[i]);
    return sqrt(dist);
}


/*!
\brief Se encarga de copiar todos los elementos del vector src al vector dst
\param src Vector de entrada.
\param dst Vector copia de src.
*/
void vectorCopy(Vector * src, Vector * dst)
{
    for (int i=0;i<src->size;i++)
        dst->data[i]=src->data[i];

}


/*!
\brief Se encarga de copiar todos los elementos de una matriz src a la matriz dst.
\param src Matriz de entrada.
\param dst Matriz copia de src.
*/
void Matrix_copy(Matrix *src, Matrix * dst)
{
    for (int i=0;i<src->size;i++)
        dst->data[i]=src->data[i];
}


/*!
\brief Se encarga de inicializar las entradas del vector src con el valor s.
\param src Vector a inicializar.
\param s Scalar.
*/
void vectorFill(Vector * src, double s)
{
    for (int i=0;i<src->size;i++)
        src->data[i]=s;
}


/*!
\brief Se encarga de inicializar las entradas de la matriz src con el valor s.
\param src Matriz a inicializar.
\param s Scalar.
*/
void matrixFill(Matrix * src, double s)
{
    for (int i=0;i<src->size;i++)
            src->data[i]=s;
}


/*!
\brief Se encarga de intercambiar elementos de un vector (v[i] = v[j], v[j] = v[i])
\param v Vector de entrada.
\param i Indice del primer elemento a intercambiar.
\param j Indice del segundo elemento a intercambiar.
*/
void vectorSwapElements(Vector *v, int i,int j)
{
}


/*!
\brief Se encarga de intercambiar todos elementos de los vectores v y u
\param v Primer vector a intercambiar.
\param u Segundo vector a intercambiar.
*/
void vectorSwap(Vector * v, Vector * u)
{
}

///Lee un sistema de ecuaciones a partir de un archivo
///el archivo debe contener el tamano del sistema, luego la matriz A y luego el vector b
void ReadEqSystem(char* filename, Matrix **A, Vector **b,Vector **x)
{
    FILE *input=fopen(filename,"r");
    if (!input){printf("File could not be opened");exit(-1);}

    int n;
    fscanf(input,"%d",&n);
    *A=RequestMatrixDoubleMem(n,n);
    *b=RequestVectorDoubleMem(n);
    *x=RequestVectorDoubleMem(n);

    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
            fscanf(input,"%lf",&(*A)->data[i*n+j]);
    for (int i=0;i<n;i++)
        fscanf(input,"%lf",&(*b)->data[i]);
}

