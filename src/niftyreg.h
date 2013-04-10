#ifndef NIFTYREG_H
#define NIFTYREG_H

typedef struct {
    nifti_image *image;
    mat44 *affine;
    int *completedIterations;
} aladin_result;

typedef struct {
    mat44 *initAffine;
    nifti_image *forwardImage;
    nifti_image *forwardControlPoints;
    nifti_image *reverseImage;
    nifti_image *reverseControlPoints;
    int *completedIterations;
} f3d_result;

extern "C"
SEXP reg_aladin_R (SEXP source, SEXP target, SEXP type, SEXP nLevels, SEXP maxIterations, SEXP useBlockPercentage, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP verbose, SEXP estimateOnly);

extern "C"
SEXP reg_f3d_R (SEXP source, SEXP target, SEXP nLevels, SEXP maxIterations, SEXP nBins, SEXP bendingEnergyWeight, SEXP jacobianWeight, SEXP inverseConsistencyWeight, SEXP finalSpacing, SEXP finalInterpolation, SEXP targetMask, SEXP sourceMask, SEXP affineComponents, SEXP initControl, SEXP symmetric, SEXP verbose, SEXP estimateOnly);

void convert_and_insert_image (nifti_image *image, SEXP list, int index);

void convert_and_insert_xform (nifti_image *image, SEXP list, int index);

void convert_and_insert_affine (mat44 *affine, SEXP list, int index);

nifti_image * s4_image_to_struct (SEXP object);

nifti_image * copy_complete_nifti_image (nifti_image *source);

mat44 * create_init_affine (nifti_image *sourceImage, nifti_image *targetImage);

nifti_image * get_deformation_field (nifti_image *targetImage, nifti_image *controlPointImage, mat44 *affineTransformation);

nifti_image * resample_image (nifti_image *sourceImage, nifti_image *targetImage, nifti_image *controlPointImage, mat44 *affineTransformation, int finalInterpolation);

aladin_result do_reg_aladin (nifti_image *sourceImage, nifti_image *targetImage, int type, int nLevels, int maxIterations, int useBlockPercentage, int finalInterpolation, nifti_image *targetMaskImage, mat44 *affineTransformation, bool verbose, bool estimateOnly);

f3d_result do_reg_f3d (nifti_image *sourceImage, nifti_image *targetImage, int nLevels, int maxIterations, int finalInterpolation, nifti_image *sourceMaskImage, nifti_image *targetMaskImage, nifti_image *controlPointImage, mat44 *affineTransformation, int nBins, float *spacing, float bendingEnergyWeight, float jacobianWeight, float inverseConsistencyWeight, bool symmetric, bool verbose, bool estimateOnly);

#endif
