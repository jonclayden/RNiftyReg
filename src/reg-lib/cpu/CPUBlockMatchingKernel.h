#ifndef CPUBLOCKMATCHINGKERNEL_H
#define CPUBLOCKMATCHINGKERNEL_H

#include "BlockMatchingKernel.h"
#include "_reg_blockMatching.h"
#ifdef RNIFTYREG
#include "RNifti.h"
#else
#include "nifti1_io.h"
#endif
#include "Content.h"

class CPUBlockMatchingKernel : public BlockMatchingKernel {
public:

    CPUBlockMatchingKernel(Content *con, std::string name) : BlockMatchingKernel(name) {
        reference = con->getCurrentReference();
        warped = con->getCurrentWarped();
        params = con->getBlockMatchingParams();
        mask = con->getCurrentReferenceMask();
    }

    void calculate();

    nifti_image *reference;
    nifti_image *warped;
    _reg_blockMatchingParam* params;
    int *mask;

};

#endif // CPUBLOCKMATCHINGKERNEL_H
