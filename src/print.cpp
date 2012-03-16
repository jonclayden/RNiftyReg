#include "print.h"

void rniftyreg_fprintf (FILE *stream, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    
    if (stream == stdout)
        Rvprintf(format, args);
    else if (stream == stderr)
        REvprintf(format, args);
    else
        vfprintf(stream, format, args);
    
    va_end(args);
}
