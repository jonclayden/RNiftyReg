#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_

nifti_image * retrieveImageFromNiftiS4 (const Rcpp::RObject &object, const bool copyData);

nifti_image * retrieveImageFromArray (const Rcpp::RObject &object);

nifti_image * retrieveImage (const SEXP _image, const bool readData);

#endif
