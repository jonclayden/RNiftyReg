#ifndef CPUOPTIMISEKERNEL_H
#define CPUOPTIMISEKERNEL_H

#include "OptimiseKernel.h"
#include "_reg_blockMatching.h"
#ifdef RNIFTYREG
#include "RNifti.h"
#else
#include "nifti1_io.h"
#endif
#include "Content.h"

class CPUOptimiseKernel : public OptimiseKernel {
public:
    CPUOptimiseKernel(Content *con, std::string name);

    _reg_blockMatchingParam *blockMatchingParams;
    mat44 *transformationMatrix;

    void calculate(bool affine, bool ils, bool svd=0);

};

#endif // CPUOPTIMISEKERNEL_H
