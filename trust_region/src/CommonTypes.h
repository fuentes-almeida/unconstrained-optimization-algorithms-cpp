#ifndef COMMONTYPES_H_
#define COMMONTYPES_H_

#define dim 6
#define imgsize 254
#define PI 3.14159265358979323846

/*!
\brief Modulo exclusivo para declarar tipo de
datos comunes entre los diferentes metodos de optimizacion.
*/


/*!
Estructura encargada de emular los elementos de un
vector n-dimensional
*/
struct Vector {
    int size;
    double *data;
};


/*!
Estructura encargada de emular los elementos de una
matriz rectangular.
*/
struct Matrix {
    int rows,
	cols,
	size;

    double *data;
};


/*!
Estructura encargada de emular los elementos de un
cubo de distintos lados.
*/
struct Cube {
    int rows,
	cols,
	depth,
	size;

    double *data;
};


/*!
\brief apuntador a funcion cuyo prototipo realiza la evaluacion de una funcion escalar R^n--->R.
\param x Vector a evaluar en la funcion.
\param pointer Apuntador Generico (opcional) cuya estructura, tipo o clase es definido por el usuario.
*/
typedef double (*scalarFunction)(struct Vector * x, void * pointer);


/*!
\brief apuntador a funcion cuyo prototipo realiza la evaluacion de una funcion vectorial R^n--->R^m.
\param x Vector a evaluar en la funcion (elemento del dominio).
\param y Evaluacion del vector (imagen de la funcion).
\param pointer Apuntador Generico (opcional) cuya estructura, tipo o clase es definido por el usuario.
*/
typedef void (*vectorialFunction)(struct Vector * x, struct Vector * y, void * pointer);


/*!
\brief apuntador a funcion cuyo prototipo realiza la evaluacion de una funcion vectorial R^n--->R^(mxn).
\param x Vector a evaluar en la funcion (elemento del dominio).
\param y Evaluacion de la matriz (imagen de la funcion).
\param pointer Apuntador Generico (opcional) cuya estructura, tipo o clase es definido por el usuario.
*/
typedef void (*matricialFunction)(struct Vector * x, struct Matrix * y, void * pointer);

typedef double (*funcion)(Matrix*,Matrix*,Matrix*,Matrix*,Vector*);


#endif
