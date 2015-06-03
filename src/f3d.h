#ifndef _F3D_H_
#define _F3D_H_

#include "nifti1_io.h"
#include "AffineMatrix.h"

struct F3dResult {
    nifti_image *forwardImage;
    nifti_image *forwardControlPoints;
    nifti_image *reverseImage;
    nifti_image *reverseControlPoints;
    std::vector<int> iterations;
};

F3dResult regF3d (nifti_image *sourceImage, nifti_image *targetImage, const int nLevels, const int maxIterations, const int interpolation, nifti_image *sourceMaskImage, nifti_image *targetMaskImage, nifti_image *controlPointImage, AffineMatrix *initAffine, const int nBins, const std::vector<float> &spacing, const float bendingEnergyWeight, const float jacobianWeight, const float inverseConsistencyWeight, const bool symmetric, const bool verbose, const bool estimateOnly);

#endif
