#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_

#include "nifti1_io.h"

nifti_image * retrieveImageFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);

nifti_image * retrieveImageFromArray (const Rcpp::RObject &object);

nifti_image * retrieveImage (const SEXP _image, const bool readData = true);

nifti_image * copyCompleteImage (const nifti_image *source);

#endif
