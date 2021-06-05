import numpy as np
import sys
import time
alpha,gamma,rho,sigma,epsilon = 1,2,0.5,0.5,1e-4
f = open('output.dat','w')

def objetiveFunction(x):
    fx=0
    for i in range(0, x.size-1):
        fx+=100*(x[i+1]-x[i]*x[i])**2+(1-x[i])**2

    return fx

def NelderMead(x,dim):

    norm,k = 1000,0
    #se ordenan los puntos de acuerdo
    #al valor de la funcion objetivo
    x=x[np.argsort(x[:,dim])]

    #centroide de los n mejores puntos
    x0=x[0:dim].mean(axis=0)

    while (norm>epsilon):
        #reflexion
        xr=x0+alpha*(x0-x[dim])

        #si el punto reflejado es mejor que el segundo peor
        # pero no mejor que el mejor, se reemplaza al peor
        xr[dim]=objetiveFunction(xr[0:dim])
        if x[0][dim] <= xr[dim] and xr[dim] < x[dim-1][dim]:
            x[dim]=np.copy(xr)

        #si el punto reflejado es el mejor
        #expansion
        elif xr[dim] < x[0][dim]:
            xe=x0+gamma*(xr-x0)
            xe[dim]=objetiveFunction(xe[0:dim])
            #si el punto expandido es mejor que el reflejado
            #se reemplaza el peor por el expandido
            if xe[dim] < xr[dim]:
                x[dim]=np.copy(xe)
            else:
                x[dim]=np.copy(xr)

        #si el punto reflejado no es mejor que el segundo peor
        #contraccion
        else:
            xc=x0+rho*(x[dim]-x0)
            xc[dim]=objetiveFunction(xc[0:dim])
            #si el punto contraido es mejor que el peor
            #se reemplaza el peor con el contraido
            if xc[dim] < x[dim][dim]:
                x[dim]=np.copy(xc)
            #si no, se reemplazan todos los puntos menos el mejor
            #reduccion
            else:
                for i in range(1,dim+1):
                    x[i]=x[0]+sigma*(x[i]-x[0])

        #se ordenan los puntos de acuerdo
        #al valor de la funcion objetivo
        x=x[np.argsort(x[:,dim])]

        #centroide de los n mejores puntos
        x0=x[0:dim].mean(axis=0)

        #se calcula la norma de la diferencia entre
        #el peor punto y el promedio de los restantes
        norm=np.linalg.norm(x[dim][0:dim]-x0[0:dim])
        k+=1

    f.write(str(k)+"\t"+str(x[0])+"\n")
    print(k,x[0])

#main funtion
#run: python Newton.py 10
if __name__ == '__main__':

    start = time.time()
    dim = int(sys.argv[1])

    #double np.array
    x = np.zeros((dim+1,dim+1))

    for i in range(0,30):
        #Se propone un punto inicial x0
        for i in range(0,dim):
            x[0][i] = np.random.uniform(-5.12,5.12,1)

        #se forma el simplejo inicial con los demas puntos
        for i in range(1,dim+1):
            x[i] = np.copy(x[0])
            x[i][i-1] += np.random.uniform(-1,1,1)

        #se evalua la funcion en los dim + 1 puntos
        for i in range(0,dim+1):
            x[i][dim]=objetiveFunction(x[i][0:dim])

        NelderMead(x,dim)
    f.write("--- %s seconds ---" % (time.time() - start))
    print("--- %s seconds ---" % (time.time() - start))
    f.close()
