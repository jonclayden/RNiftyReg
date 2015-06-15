#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_

#include "nifti1_io.h"

// Thin wrapper around a C-style nifti_image struct that allows C++-style destruction
class NiftiImage
{
protected:
    nifti_image *image;
    
    void copy (const NiftiImage &source)
    {
        nifti_image *sourceStruct = source;
        size_t dataSize = nifti_get_volsize(sourceStruct);
        image = nifti_copy_nim_info(sourceStruct);
        image->data = calloc(1, dataSize);
        memcpy(image->data, sourceStruct->data, dataSize);
    }
    
public:
    NiftiImage ()
        : image(NULL) {}
    
    NiftiImage (const NiftiImage &source)
    {
        copy(source);
    }
    
    NiftiImage (nifti_image * const image)
        : image(image) {}
    
    ~NiftiImage ()
    {
        nifti_image_free(image);
    }
    
    operator nifti_image* () const { return image; }
    
    NiftiImage & operator= (const NiftiImage &source)
    {
        copy(source);
        return *this;
    }
    
    bool isNull () const { return (image == NULL); }
};

NiftiImage retrieveImageFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);

NiftiImage retrieveImageFromArray (const Rcpp::RObject &object);

NiftiImage retrieveImage (const SEXP _image, const bool readData = true);

Rcpp::RObject imageToArray (nifti_image *source);

#endif
