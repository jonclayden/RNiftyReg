#include <Rmath.h>

#include "index.h"

void matrix_to_vector_loc (const int *loc, const int *dim, const int n_dims, size_t *vector_loc)
{
    int i, j;
    size_t temp;
    
    *vector_loc = loc[0];
    
    for (i=1; i<n_dims; i++)
    {
        temp = loc[i];
        for (j=0; j<i; j++)
            temp *= dim[j];
        *vector_loc += temp;
    }
}

void vector_to_matrix_loc (const size_t loc, const int *dim, const int n_dims, int *matrix_loc)
{
    int i, j;
    size_t temp;
    
    matrix_loc[0] = loc % dim[0];
    
    for (i=1; i<n_dims; i++)
    {
        temp = 1;
        for (j=0; j<i; j++)
            temp *= dim[j];
        matrix_loc[i] = (loc / temp) % dim[i];
    }
}

unsigned char index_uchar_array (const unsigned char *array, const int *loc, const int *dim, const int n_dims)
{
    size_t vector_loc;
    matrix_to_vector_loc(loc, dim, n_dims, &vector_loc);
    return array[vector_loc];
}

int index_int_array (const int *array, const int *loc, const int *dim, const int n_dims)
{
    size_t vector_loc;
    matrix_to_vector_loc(loc, dim, n_dims, &vector_loc);
    return array[vector_loc];
}

float index_float_array (const float *array, const int *loc, const int *dim, const int n_dims)
{
    size_t vector_loc;
    matrix_to_vector_loc(loc, dim, n_dims, &vector_loc);
    return array[vector_loc];
}

double index_double_array (const double *array, const int *loc, const int *dim, const int n_dims)
{
    size_t vector_loc;
    matrix_to_vector_loc(loc, dim, n_dims, &vector_loc);
    return array[vector_loc];
}

int loc_in_bounds (const int *loc, const int *dim, const int n_dims)
{
    int in_bounds = 1;
    
    for (int i=0; i<n_dims; i++)
    {
        if (loc[i] < 0 || loc[i] > (dim[i]-1))
        {
            in_bounds = 0;
            break;
        }
    }
    
    return in_bounds;
}

double inner_product (const double *a, const double *b, const int len)
{
    double result = 0.0;
    
    for (int i=0; i<len; i++)
        result += a[i] * b[i];
    
    return (result);
}

double vector_length (const double *a, const int len)
{
    double length = 0.0;
    
    for (int i=0; i<len; i++)
        length += R_pow_di(a[i], 2);
    
    return (sqrt(length));
}
