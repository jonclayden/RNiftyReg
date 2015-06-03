#ifndef _ALADIN_H_
#define _ALADIN_H_

#include "nifti1_io.h"
#include "AffineMatrix.h"

enum LinearTransformScope { RigidScope, AffineScope };

struct AladinResult
{
    nifti_image *image;
    AffineMatrix *affine;
    std::vector<int> iterations;
};

AladinResult regAladin (nifti_image *sourceImage, nifti_image *targetImage, const LinearTransformScope scope, const bool symmetric, const int nLevels, const int maxIterations, const int useBlockPercentage, const int interpolation, nifti_image *sourceMaskImage, nifti_image *targetMaskImage, AffineMatrix *initAffine, const bool verbose, const bool estimateOnly);

#endif
