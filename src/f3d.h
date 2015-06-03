#ifndef _F3D_H_
#define _F3D_H_

#include "nifti1_io.h"
#include "AffineMatrix.h"

struct f3dResult {
    // AffineMatrix *initAffine;
    nifti_image *forwardImage;
    nifti_image *forwardControlPoints;
    nifti_image *reverseImage;
    nifti_image *reverseControlPoints;
    std::vector<int> iterations;
};

#endif
