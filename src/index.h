#ifndef _INDEX_H_
#define _INDEX_H_

#include <stddef.h>

void matrix_to_vector_loc (const int *loc, const int *dim, const int n_dims, size_t *vector_loc);

void vector_to_matrix_loc (const size_t loc, const int *dim, const int n_dims, int *matrix_loc);

unsigned char index_uchar_array (const unsigned char *array, const int *loc, const int *dim, const int n_dims);

int index_int_array (const int *array, const int *loc, const int *dim, const int n_dims);

float index_float_array (const float *array, const int *loc, const int *dim, const int n_dims);

double index_double_array (const double *array, const int *loc, const int *dim, const int n_dims);

int loc_in_bounds (const int *loc, const int *dim, const int n_dims);

double inner_product (const double *a, const double *b, const int len);

double vector_length (const double *a, const int len);

#endif
