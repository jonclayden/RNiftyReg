#ifndef _F3D_H_
#define _F3D_H_

#include "RNifti.h"
#include "AffineMatrix.h"

struct F3dResult
{
    NiftiImage image;
    NiftiImage forwardTransform;
    NiftiImage reverseTransform;
    std::vector<int> iterations;
};

F3dResult regF3d (const NiftiImage &sourceImage, const NiftiImage &targetImage, const int nLevels, const int maxIterations, const int interpolation, const NiftiImage &sourceMaskImage, const NiftiImage &targetMaskImage, const NiftiImage &initControlPoints, const AffineMatrix &initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float linearEnergyWeight, const float jacobianWeight, const bool symmetric, const bool verbose, const bool estimateOnly);

#endif
