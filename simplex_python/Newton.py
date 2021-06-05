import numpy as np
import sys

def objetiveFunction(x):

	z = x[0]**2 + x[1]**2

	return z


def gradientFunction(x):

	grad = np.zeros( x.shape )

	grad[0] = 2*x[0]
	grad[1] = 2*x[1]

	return grad


def hessianFunction(x):

	hessian = np.zeros( (x.shape[0], x.shape[0]) )

	hessian[0][0] = 2
	hessian[1][1] = 2

	return hessian



#Newton solver.
#funtion Objetive Function.
#gradient Gradient of Function.
#hessian Matrix hessian.
#x0 Initial value.
#iterMax Maximum of iterations
def solver(function, gradient, hessian, x0, iterMax):

	alpha = 0.5
	xk = np.copy(x0)

	for k in range(0, iterMax):

		g = -1.0*gradient(xk)
		H = hessian(xk)
		pk = np.linalg.solve(H, g)
		xk = xk + alpha*pk

		print "k: ", k, ", f(xk): " ,function(xk)

	return xk


#main funtion
#run: python Newton.py 10
if __name__ == '__main__':

	#arguments
	argc = len(sys.argv)

	#simple np.array
	x0 = np.zeros((2))

	#random init
	x0[0] = np.random.uniform(-4,4,1)
	x0[0] = np.random.uniform(-4,4,1)

	#size of array
	print "dim x: \n", x0.shape, "\n"

	#value of x0
	print "x0: \n ", x0, "\n"

	#alias or "pointers of funtion"
	function = objetiveFunction
	gradient = gradientFunction
	hessian = hessianFunction

	print "f(x): \n", function(x0), "\n"
	print "Gf(x): \n",  gradient(x0), "\n"
	print "Hf(x): \n", hessian(x0), "\n"

	#solver
	x = solver(function, gradient, hessian, x0, int(sys.argv[1]))
	print "\nsolution: \n ", x, "\n"
	
