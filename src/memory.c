#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <assert.h>

#include "memory.h"

//for allocating an array of dimensions dim1*dim2 and of bytes size si
void** calloc_2D(int dim1, int dim2, size_t si)
{
    int i;
    void **array=NULL;

    //Here, allocation of the first dimension : rows
    array = calloc(dim1,si);
    assert(array!=NULL);

    //Then allocation of the second dimensions : columns
    for (i=0; i<dim1; i++)
    {
        array[i] = calloc(dim2,si);
        assert(array[i]!=NULL);
    }

    return array;
}

//for freeing one or more dynamically allocated 2D arrays of first dimension dim1
void free_2D(int dim1, ...)
{
    int i;
    void **array=NULL;
    va_list ap;

    va_start(ap,dim1);
    array=va_arg(ap,void**);

    while(array!=NULL)
    {
        for(i=0; i<dim1; i++)
            free(array[i]);
        free(array);
        array=NULL;
        array=va_arg(ap,void**);
    }
    va_end(ap);
}

//for allocating an array of dimensions dim1*dim2*dim3 and of bytes size si
void*** calloc_3D(int dim1, int dim2, int dim3, size_t si)
{
    int i,j;
    void ***array=NULL;

    array=calloc(dim1,si);
    assert(array!=NULL);

    for(i=0; i<dim1; i++)
    {
        array[i]=calloc(dim2,si);
        assert(array[i]!=NULL);

        for(j=0; j<dim2; j++)
        {
            array[i][j]=calloc(dim3,si);
            assert(array[i]!=NULL);
        }
    }
    return array;
}

//for freeing one or more dynamically allocated 3D arrays of dimensions dim1,dim2
void free_3D(int dim1, int dim2, ...)
{
    int i,j;
    void ***array=NULL;
    va_list ap;

    va_start(ap,dim2);
    array=va_arg(ap,void***);

    while(array!=NULL)
    {
        for(i=0; i<dim1; i++)
        {
            for(j=0; j<dim2; j++)
                free(array[i][j]);
            free(array[i]);
        }
        free(array);
        array=NULL;
        array=va_arg(ap,void***);
    }
    va_end(ap);
}
