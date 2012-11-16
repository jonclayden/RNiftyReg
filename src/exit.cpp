#include <R.h>

#include "exit.h"

void rniftyreg_exit (int status)
{
    error("NiftyReg encountered a fatal error (status %d)", status);
}
