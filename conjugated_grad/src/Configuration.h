#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "CommonTypes.h"

/*!
\brief Modulo exclusivo para contener variables
globales de configuracion que se emplean en los
diferentes metodos de optimizacion y que se inicializan
dinamicamente en cada problema a resolver.

Cada variable al agregarse debe tener la siguiente notacion g_nombreVariable
seguida de un comentario breve sobre el uso de la variable.
*/
double c1,c2;

struct Config
{
	int maxIter; //Maximo de iteraciones a ejecutar en cada metodo
	double eps; //Epsilon del sistema
};

extern Config g_config;

#endif
