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
            
            size_t blockSize = 1;
            for (int i=1; i<dimension; i++)
                blockSize *= image->dim[i];
            
            if (blockSize != source->nvox)
                throw std::runtime_error("New data does not have the same size as the target block");
            
            blockSize *= image->nbyper;
            memcpy(static_cast<char*>(image->data) + blockSize*index, source->data, blockSize);
            return *this;
        }
    };
    
protected:
    nifti_image *image;
    bool persistent;
    
    void copy (nifti_image * const source);
    void copy (const NiftiImage &source);
    void copy (const Block &source);
    
    void initFromNiftiS4 (const Rcpp::RObject &object, const bool copyData = true);
    void initFromArray (const Rcpp::RObject &object);
    
public:
    NiftiImage ()
        : image(NULL), persistent(false) {}
    
    NiftiImage (const NiftiImage &source)
        : persistent(false)
    {
        this->copy(source);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p\n", this->image);
#endif
    }
    
    NiftiImage (nifti_image * const image, const bool copy = false)
        : image(image), persistent(false)
    {
        if (copy)
            this->copy(image);
#ifndef NDEBUG
        Rprintf("Creating NiftiImage with pointer %p\n", this->image);
#endif
    }
    
    NiftiImage (const SEXP object, const bool readData = true);
    
    ~NiftiImage ()
    {
        if (!persistent)
        {
#ifndef NDEBUG
            Rprintf("Freeing NiftiImage with pointer %p\n", this->image);
#endif
            nifti_image_free(image);
        }
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
    
    void setPersistence (const bool persistent)
    {
        this->persistent = persistent;
#ifndef NDEBUG
        if (persistent)
            Rprintf("Setting NiftiImage with pointer %p to be persistent\n", this->image);
#endif
    }
    
    bool isNull () const { return (image == NULL); }
    int nDims () const
    {
        if (image == NULL)
            return 0;
        else
            return image->ndim;
    }
    
    // Note that this function differs from its R equivalent in only dropping unitary dimensions after the last nonunitary one
    NiftiImage & drop ()
    {
        int ndim = image->ndim;
        while (image->dim[ndim] < 2)
            ndim--;
        image->dim[0] = image->ndim = ndim;
        
        return *this;
    }
    
    void update (const SEXP array);
    
    mat44 xform (const bool preferQuaternion = true) const;
    
    const Block slice (const int i) const { return Block(*this, 3, i); }
    const Block volume (const int i) const { return Block(*this, 4, i); }
    
    Block slice (const int i) { return Block(*this, 3, i); }
    Block volume (const int i) { return Block(*this, 4, i); }
    
    Rcpp::RObject toArray () const;
    Rcpp::RObject toPointer (const std::string label) const;
    Rcpp::RObject headerToList () const;
};

NiftiImage allocateMultiregResult (const NiftiImage &source, const NiftiImage &target, const bool forceDouble);

#endif
