#ifndef _NIFTI_IMAGE_H_
#define _NIFTI_IMAGE_H_

#include "nifti1_io.h"

// Thin wrapper around a C-style nifti_image struct that allows C++-style destruction
class NiftiImage
{
public:
    struct Block
    {
        const NiftiImage &image;
        const int dimension;
        const int index;
        
        Block (const NiftiImage &image, const int dimension, const int index)
            : image(image), dimension(dimension), index(index)
        {
            if (dimension != image->ndim)
                throw std::runtime_error("Blocks must be along the last dimension in the image");
        }
        
        Block & operator= (const NiftiImage &source)
        {
            if (source->datatype != image->datatype)
                throw std::runtime_error("New data does not have the same datatype as the target block");
            
            size_t blockSize = 0;
            for (int i=1; i<dimension; i++)
                blockSize += source->dim[i];
            memcpy(static_cast<char*>(image->data) + blockSize*index, source->data, blockSize);
            return *this;
        }
    };
    
protected:
    nifti_image *image;
    bool persistent;
    
    void copy (nifti_image * const source)
    {
        if (source != NULL)
        {
            size_t dataSize = nifti_get_volsize(source);
            image = nifti_copy_nim_info(source);
            image->data = calloc(1, dataSize);
            memcpy(image->data, source->data, dataSize);
        }
    }
    
    void copy (const NiftiImage &source)
    {
        nifti_image *sourceStruct = source;
        copy(sourceStruct);
    }
    
    void copy (const Block &source)
    {
        nifti_image *sourceStruct = source.image;
        if (sourceStruct != NULL)
        {
            image = nifti_copy_nim_info(sourceStruct);
            image->dim[0] = source.image->dim[0] - 1;
            image->dim[source.dimension] = 1;
            image->pixdim[source.dimension] = 1.0;
            nifti_update_dims_from_array(image);
            
            size_t blockSize = 0;
            for (int i=1; i<source.dimension; i++)
                blockSize += source.image->dim[i];
            image->data = calloc(1, blockSize);
            memcpy(image->data, static_cast<char*>(source.image->data) + blockSize*source.index, blockSize);
        }
    }
    
public:
    NiftiImage ()
        : image(NULL), persistent(false) {}
    
    NiftiImage (const NiftiImage &source)
        : persistent(false)
    {
        this->copy(source);
    }
    
    NiftiImage (nifti_image * const image, const bool copy = false)
        : image(image), persistent(false)
    {
        if (copy)
            this->copy(image);
    }
    
    ~NiftiImage ()
    {
        if (!persistent)
            nifti_image_free(image);
    }
    
    operator nifti_image* () const { return image; }
    
    nifti_image * operator-> () const { return image; }
    
    NiftiImage & operator= (const NiftiImage &source)
    {
        copy(source);
        return *this;
    }
    
    NiftiImage & operator= (const Block &source)
    {
        copy(source);
        return *this;
    }
    
    void setPersistence (const bool persistent) { this->persistent = persistent; }
    
    bool isNull () const { return (image == NULL); }
    int nDims () const
    {
        if (image == NULL)
            return 0;
        else
            return image->ndim;
    }
    
    const Block slice (const int i) const { return Block(*this, 3, i); }
    const Block volume (const int i) const { return Block(*this, 4, i); }
    
    Block slice (const int i) { return Block(*this, 3, i); }
    Block volume (const int i) { return Block(*this, 4, i); }
};

NiftiImage allocateMultiregResult (const NiftiImage &source, const NiftiImage &target, const bool forceDouble);

NiftiImage retrieveImageFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);

NiftiImage retrieveImageFromArray (const Rcpp::RObject &object);

NiftiImage retrieveImage (const SEXP _image, const bool readData = true);

Rcpp::RObject imageToArray (nifti_image *source);

Rcpp::RObject imageToPointer (nifti_image *source, const std::string label);

#endif
