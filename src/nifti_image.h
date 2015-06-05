#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_

#include "nifti1_io.h"

// Thin wrapper around a C-style nifti_image struct that allows C++-style destruction
class NiftiImage
{
protected:
    nifti_image *image;
    
public:
    NiftiImage ()
        : image(NULL) {}
    
    NiftiImage (nifti_image *image)
        : image(image) {}
    
    ~NiftiImage ()
    {
        nifti_image_free(image);
    }
    
    operator nifti_image* () const
    {
        return image;
    }
};

nifti_image * retrieveImageFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);

nifti_image * retrieveImageFromArray (const Rcpp::RObject &object);

nifti_image * retrieveImage (const SEXP _image, const bool readData = true);

nifti_image * copyCompleteImage (const nifti_image *source);

Rcpp::RObject imageToArray (nifti_image *source);

#endif
