#ifndef GAUSSNEWTON_H_INCLUDED
#define GAUSSNEWTON_H_INCLUDED

void GaussNewtonMethod(Vector **Jacobian, Matrix *JtJ, Vector *Rx, Vector *JR, Vector *pk, int dim)
{
    matrixFill(JtJ,0);
    vectorFill(JR,0);

    for (int i=0;i<dim;i++)
        for (int j=i;j<dim;j++)
        {
            JtJ->data[i*dim+j]=vectorInnerProduct(Jacobian[i],Jacobian[j]);
            if (i!=j) JtJ->data[j*dim+i]=JtJ->data[i*dim+j];
        }

    for (int k=0;k<dim;k++)
        JR->data[k]=-vectorInnerProduct(Jacobian[k],Rx);

    //ConjugatedGradient(JtJ,JR,pk);
    //GaussianElimination(JtJ,JR,pk,dim);
    CholeskySolve(JtJ,JR,pk);
    //SolveCrouthSystem(JtJ,JR,pk,dim);
}


#endif // GAUSSNEWTON_H_INCLUDED

