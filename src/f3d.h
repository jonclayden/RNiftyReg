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
};

F3dResult regF3d (RNifti::NiftiImage &sourceImage, RNifti::NiftiImage &targetImage, const int nLevels, const int maxIterations, const int interpolation, RNifti::NiftiImage &sourceMaskImage, RNifti::NiftiImage &targetMaskImage, RNifti::NiftiImage &initControlPoints, const AffineMatrix &initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float linearEnergyWeight, const float jacobianWeight, const bool symmetric, const bool verbose, const bool estimateOnly);

#endif
