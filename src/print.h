#ifndef PRINT_H
#define PRINT_H

#define R_USE_C99_IN_CXX
#include <R_ext/Print.h>

#include <stdio.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C"
#endif
void rniftyreg_fprintf (FILE *stream, const char *format, ...);

#ifdef __cplusplus
extern "C"
#endif
int rniftyreg_fputs (const char *str, FILE *stream);

#endif
