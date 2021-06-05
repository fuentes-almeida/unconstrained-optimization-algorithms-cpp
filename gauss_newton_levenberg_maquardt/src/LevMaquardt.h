#ifndef LEVMAQUARDT_H_INCLUDED
#define LEVMAQUARDT_H_INCLUDED

void LevenbergMaquardtMethod(Vector **Jacobian, Matrix *JtJ, Vector *Rx, Vector *JR, Vector *pk, int dim)
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

    Vector *Aux1=RequestVectorDoubleMem(dim); vectorFill(Aux1,1);
    Matrix *Aux2=RequestMatrixDoubleMem(dim,dim); Matrix_copy(JtJ,Aux2);

    double lambda_min=InversePowerMethod(Aux2,Aux1,dim,1e-6,1000);

    FreeVectorMem(Aux1);
    FreeMatrixMem(Aux2);

    for (int i=0;i<dim;i++)
        JtJ->data[i*dim+i]+=lambda_min;

    ConjugatedGradient(JtJ,JR,pk);
    //GaussianElimination(JtJ,JR,pk,dim);
    //CholeskySolve(AlphaJtJ,AlphaJR,AlphaPk);
    //SolveCrouthSystem(JtJ,JR,pk,dim);
}


#endif // LEVMAQUARDT_H_INCLUDED
