#ifndef _ALADIN_H_
#define _ALADIN_H_

#include "nifti_image.h"
#include "AffineMatrix.h"

enum LinearTransformScope { RigidScope, AffineScope };

struct AladinResult
{
    NiftiImage image;
    AffineMatrix affine;
    std::vector<int> iterations;
};

AladinResult regAladin (const NiftiImage &sourceImage, const NiftiImage &targetImage, const LinearTransformScope scope, const bool symmetric, const int nLevels, const int maxIterations, const int useBlockPercentage, const int interpolation, const NiftiImage &sourceMaskImage, const NiftiImage &targetMaskImage, const AffineMatrix &initAffine, const bool verbose, const bool estimateOnly);

#endif
