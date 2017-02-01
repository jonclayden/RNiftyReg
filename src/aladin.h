#ifndef _ALADIN_H_
#define _ALADIN_H_

#include "RNifti.h"
#include "AffineMatrix.h"

enum LinearTransformScope { RigidScope, AffineScope };

struct AladinResult
{
    RNifti::NiftiImage image;
    AffineMatrix forwardTransform;
    AffineMatrix reverseTransform;
    std::vector<int> iterations;
};

AladinResult regAladin (RNifti::NiftiImage &sourceImage, RNifti::NiftiImage &targetImage, const LinearTransformScope scope, const bool symmetric, const int nLevels, const int maxIterations, const int useBlockPercentage, const int interpolation, RNifti::NiftiImage &sourceMaskImage, RNifti::NiftiImage &targetMaskImage, const AffineMatrix &initAffine, const bool verbose, const bool estimateOnly);

#endif
