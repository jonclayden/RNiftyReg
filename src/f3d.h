#ifndef _F3D_H_
#define _F3D_H_

#include "RNifti.h"
#include "AffineMatrix.h"

struct F3dResult
{
    RNifti::NiftiImage image;
    RNifti::NiftiImage forwardTransform;
    RNifti::NiftiImage reverseTransform;
    std::vector<int> iterations;
    RNifti::NiftiImage source;
    RNifti::NiftiImage target;
};

template <typename PrecisionType>
F3dResult regF3d (const RNifti::NiftiImage &sourceImage, const RNifti::NiftiImage &targetImage, const int nLevels, const int maxIterations, const int interpolation, const RNifti::NiftiImage &sourceMaskImage, const RNifti::NiftiImage &targetMaskImage, const RNifti::NiftiImage &initControlPoints, const AffineMatrix &initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float linearEnergyWeight, const float jacobianWeight, const bool symmetric, const bool verbose, const bool estimateOnly);

#endif
