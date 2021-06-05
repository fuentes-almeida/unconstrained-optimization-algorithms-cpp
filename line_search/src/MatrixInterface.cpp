#include <stdlib.h>
#include <math.h>
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
}

/*!
\brief Se encarga de realizar la operaci\'on vectorial c = a + b
\param a Vector de entrada
\param b Vector de entrada
\param c Vector de resultado
*/
void vectorAdd(Vector *a, Vector *b, Vector *c)
{
    for (int i=0;i<dim;i++)
        c->data[i]=a->data[i]+b->data[i];
}


/*!
\brief Se encarga de realizar la operaci\'on vectorial b = a*s
\param a Vector de entrada
\param b Vector de resultado
\param s Scalar
*/
void vectorScalarProduct(Vector *a, Vector * b, double s)
{
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
	return 0;
}


/*!
\brief Se encarga de calcular la norma  del vector, norm = ||a||
\param a Vector de entrada.
\return n Norma del vector.
*/
double vectorNorm(Vector *a)
{
	return 0;
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
    double result=0;
    for (int i=0;i<dim;i++)
        result+=a->data[i]*b->data[i];

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
}


/*!
\brief Se encarga de realizar la operaci\'on matricial C = A - B
\param A Matriz de entrada
\param B Matriz de entrada
\param C Matriz resultante
*/
void matrixSubstract(Matrix * A, Matrix * B, Matrix * C)
{
}


/*!
\brief Se encarga de realizar la operaci\'on matricial B = A*s
\param A Matriz de entrada
\param B Matriz resultante
\param s Scalar
*/
void matrixScalarProduct(Matrix * A, Matrix * B, double s)
{

}


/*!
\brief Se encarga de realizar la operaci\'on sobre la diagonal de A,  B(i,i) = A(i,i) + s
\param A Matriz de entrada
\param B Matriz resultante
\param s Scalar
*/
void matrixAddToDiagonal(Matrix * A, Matrix * B, double s)
{
}


/*!
\brief Se encarga de realizar la operaci\'on matricial B = A^T
\param a Vector de entrada
\param b Vector de resultado
\param s Scalar
*/
void matrixTranspose(Matrix * A, Matrix * B)
{
}



/*!
\brief Se encarga de realizar la operaci\'on matricial C = A*B
\param A Matriz de entrada
\param B Matriz de entrada
\param C Matriz resultante
*/
void matrixProduct(Matrix * A, Matrix * B, Matrix * C)
{
}


/*!
\brief se encarga de encontrar el m\'inimo elemento de la matriz A.
\param A Matriz a calcular el m\'inimo.
\return b M\'inimo elemento de la matriz.
*/
double matrixMin(Matrix *A)
{
	return 0;
}


/*!
\brief se encarga de encontrar el m\'aximo elemento de la matriz A.
\param A Matriz a calcular el m\'aximo.
\return b M\'aximo elemento de la matriz.
*/
double matrixMax(Matrix *A)
{
	return 0;
}



/*!
\brief se encarga de calcular el producto entre una matriz y un vector, b = A * x.
\param A Matriz de entrada.
\param x Vector de incognitas.
\param b Vector resultante de la multiplicaci\'on.
*/

void matrixVectorProduct(Matrix * A, Vector * x, Vector * b)
{
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
	return 0;
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

void Gradient(funcion fx,Vector *x,Vector *v)
{
    for (int i=0;i<dim;i++)
    v->data[i]=PartialDeriv(x,i,STEP,fx);
}

void Hessian(funcion fx,Vector *x,Matrix *matrix)
{
    for (int i=0;i<dim;i++)
        for (int j=0;j<dim;j++)
            matrix->data[i*dim+j]=MixPartialDeriv(x,i,j,STEP,STEP,fx);
}


/***********************************************************Factorizaciones y solvers para sistemas de ecuaciones Ax = b************************************/


/*!\brief Obtiene la factorizaci\'on A = LU de la matriz de entrada A, usando el m\'etodo de Crout.
\param A Matriz a descomponer.
\param L Matriz a inferior.
\param U Matriz a superior.
*/
void croutDecomposition(Matrix * A, Matrix * L, Matrix * U)
{
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
    int n=dim;
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
    int n=dim;

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
    int n=dim;

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
    int n=dim;
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


/***********************************************************************************************Funciones de utileria********************/

double Rosenbrock(Vector *x)
{
    double fx=0;
    for (int i=0;i<dim-1;i++)
        fx+=100*(x->data[i+1]-x->data[i]*x->data[i])*(x->data[i+1]-x->data[i]*x->data[i])+(1-x->data[i])*(1-x->data[i]);
    return fx;
}

double Euclidean(Vector *v1,Vector *v2)
{
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
}


/*!
\brief Se encarga de copiar todos los elementos de una matriz src a la matriz dst.
\param src Matriz de entrada.
\param dst Matriz copia de src.
*/
void Matrix_copy(Matrix src, Matrix * dst)
{
}


/*!
\brief Se encarga de inicializar las entradas del vector src con el valor s.
\param src Vector a inicializar.
\param s Scalar.
*/
void vectorFill(Vector * src, double s)
{
}


/*!
\brief Se encarga de inicializar las entradas de la matriz src con el valor s.
\param src Matriz a inicializar.
\param s Scalar.
*/
void matrixFill(Matrix * src, double s)
{
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



