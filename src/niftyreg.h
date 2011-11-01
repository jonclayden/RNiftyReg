#ifndef NIFTYREG_H
#define NIFTYREG_H

typedef struct {
    nifti_image *image;
    mat44 *affine;
    int *completedIterations;
} aladin_result;

extern "C"
SEXP reg_aladin (SEXP source, SEXP target, SEXP type, SEXP finalPrecision, SEXP nLevels, SEXP maxIterations, SEXP useBlockPercentage, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP verbose);

nifti_image * s4_image_to_struct (SEXP object);

bool reg_test_convergence (mat44 *updateMatrix);

nifti_image * copy_complete_nifti_image (nifti_image *source);

nifti_image * create_position_field (nifti_image *templateImage, bool twoDimRegistration);

aladin_result do_reg_aladin (nifti_image *sourceImage, nifti_image *targetImage, int type, int finalPrecision, int nLevels, int maxIterations, int useBlockPercentage, int finalInterpolation, nifti_image *targetMaskImage, mat44 *affineTransformation, bool verbose);

#endif
