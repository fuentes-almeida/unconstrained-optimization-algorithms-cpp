#include "PGMFunctions.h"
#include "MatrixInterface.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*!
\brief Funcion encargada de quitar los espacios en blanco
\param file Apuntador de archivo del documento abierto.
*/
void skipLine(FILE *file) {
    char c;
    do {
        c = fgetc(file);
    }
    while (c != '\n');
}


/*!
\brief Read a PGM file.
\param filename Nombre de la imagen a abrir.
\param img Matriz que contrenda los valores de la imagen.
*/
int readPGM(const char *filename, struct Matrix  * img)
{
    FILE *file = fopen(filename,"r");

    if (file==NULL)
    {
        printf("Could not open file: %s\n",filename);
	return -1;
    }

    skipLine(file);

    char c = '#';

    while (c == '#')
    {
        c = fgetc(file);
        ungetc(c, file);
        if (c == '#')  skipLine(file);
    }

    double maxValue;

    fscanf(file,"%d",&img->cols);
    fscanf(file,"%d",&img->rows);
    fscanf(file,"%lf",&maxValue);

    img->size=img->cols*img->rows;

    img->data = (double *)calloc( img->rows*img->cols, sizeof(double));

    for(int i=0; i< img->rows; i++)
        for(int j=0; j< img->cols; j++)
            fscanf(file,"%lf", &img->data[i*img->cols + j]);
    fclose(file);

    return 0;
}

/*!
\brief Write a matrix to a PGM file.
\param filename Nombre de la imagen a guardar.
\param img Matriz que contiene los valores de la imagen.
*/
void writePGM(const char *filename, struct Matrix * img)
{
    FILE *file = fopen(filename,"w");

    fprintf(file,"P2\n#Created by image.c\n%d %d\n%d\n", img->cols,
	img->rows,255);

    for(int i=0; i < img->rows; i++)
    {
        for(int j=0; j < img->cols; j++) {
            fprintf(file,"%d ",(int)img->data[i*img->cols + j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

void ImageGradient(Matrix *Image, Matrix **ImgGrad)
{
    for (int i=0;i<Image->rows;i++)
        for (int j=0;j<Image->cols;j++)
        {
            if (j==0)
                ImgGrad[0]->data[i*Image->cols+j]=0.5*(Image->data[i*Image->cols+j+1]-Image->data[i*Image->cols+Image->cols-1]);
            else if (j==Image->cols-1)
                ImgGrad[0]->data[i*Image->cols+j]=0.5*(Image->data[i*Image->cols]-Image->data[i*Image->cols+j-1]);
            else
                ImgGrad[0]->data[i*Image->cols+j]=0.5*(Image->data[i*Image->cols+j+1]-Image->data[i*Image->cols+j-1]);
        }

    for (int i=0;i<Image->rows;i++)
        for (int j=0;j<Image->cols;j++)
        {
            if (i==0)
                ImgGrad[1]->data[i*Image->cols+j]=0.5*(Image->data[(i+1)*Image->cols+j]-Image->data[(Image->cols-1)*Image->cols+j]);
            else if (i==Image->cols-1)
                ImgGrad[1]->data[i*Image->cols+j]=0.5*(Image->data[Image->cols+j]-Image->data[(i-1)*Image->cols+j]);
            else
                ImgGrad[1]->data[i*Image->cols+j]=0.5*(Image->data[(i+1)*Image->cols+j]-Image->data[(i-1)*Image->cols+j]);
        }
}

void ImageHessian(Matrix **ImgGrad, Matrix **ImgHessian)
{
    ImageGradient(ImgGrad[0],&ImgHessian[0]);
    ImageGradient(ImgGrad[1],&ImgHessian[2]);
}

void NormalizeImage(Matrix *Image, Matrix *Normalized)
{
    double maxi,mini;
    maxi=matrixMax(Image);
    mini=matrixMin(Image);

    for (int i=0;i<Image->rows;i++)
        for (int j=0;j<Image->cols;j++)
            Normalized->data[i*Image->cols+j]=255*(Image->data[i*Image->cols+j]-mini)/(maxi-mini);

}

int NumNeighbours(Matrix *Image, int row, int col)
{
    int N=Image->rows,M=Image->cols;
    int neighbours=0;

    if (row>0)
        neighbours++;
    if (row<N-1)
        neighbours++;
    if (col>0)
        neighbours++;
    if (col<M-1)
        neighbours++;

    return neighbours;
}

double PixelNeighbours(Matrix *Image, int row, int col)
{
    int N=Image->rows;
    int M=Image->cols;
    double sum=0;

    if (row>0)
        sum+=Image->data[(row-1)*M+col];
    if (row<N-1)
        sum+=Image->data[(row+1)*M+col];
    if (col>0)
        sum+=Image->data[row*M+col-1];
    if (col<M-1)
        sum+=Image->data[row*M+col+1];

    return sum;
}

double ErrorFunction(Matrix *ImageX, Matrix *ImageY,double lambda)
{
    int N=ImageX->rows,M=ImageX->cols;
    double error=0;

    for (int i=0;i<N;i++)
        for (int j=0;j<M;j++)
        {
            error+=0.5*pow(ImageX->data[i*M+j]-ImageY->data[i*M+j],2);
            if (i>0)
                error+=0.5*lambda*pow(ImageX->data[i*M+j]-ImageY->data[(i-1)*M+j],2);
            if (j>0)
                error+=0.5*lambda*pow(ImageX->data[i*M+j]-ImageY->data[i*M+(j-1)],2);
        }
    return error;
}


