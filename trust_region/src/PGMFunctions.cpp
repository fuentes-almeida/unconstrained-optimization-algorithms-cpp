#include "PGMFunctions.h"
#include "MatrixInterface.h"
#include <stdio.h>
#include <stdlib.h>

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
    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
        {
            if (j==0)
                ImgGrad[0]->data[i*imgsize+j]=0.5*(Image->data[i*imgsize+j+1]-Image->data[i*imgsize+imgsize-1]);

            else if (j==imgsize-1)
                ImgGrad[0]->data[i*imgsize+j]=0.5*(Image->data[i*imgsize]-Image->data[i*imgsize+j-1]);

            else
                ImgGrad[0]->data[i*imgsize+j]=0.5*(Image->data[i*imgsize+j+1]-Image->data[i*imgsize+j-1]);
        }

    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
        {

            if (i==0)
                ImgGrad[1]->data[i*imgsize+j]=0.5*(Image->data[(i+1)*imgsize+j]-Image->data[(imgsize-1)*imgsize+j]);
            else if (i==imgsize-1)
                ImgGrad[1]->data[i*imgsize+j]=0.5*(Image->data[imgsize+j]-Image->data[(i-1)*imgsize+j]);
            else
                ImgGrad[1]->data[i*imgsize+j]=0.5*(Image->data[(i+1)*imgsize+j]-Image->data[(i-1)*imgsize+j]);

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

    for (int i=0;i<imgsize;i++)
        for (int j=0;j<imgsize;j++)
            Normalized->data[i*imgsize+j]=255*(Image->data[i*imgsize+j]-mini)/(maxi-mini);

}




