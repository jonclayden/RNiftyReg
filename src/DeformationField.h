#ifndef _DEFORMATION_FIELD_H_
#define _DEFORMATION_FIELD_H_

#include <RcppEigen.h>

#include "config.h"
#include "nifti1_io.h"
#include "_reg_localTransformation.h"
#include "_reg_globalTransformation.h"
#include "_reg_resampling.h"
#include "NiftiImage.h"
#include "AffineMatrix.h"

class DeformationField
{
protected:
    NiftiImage deformationFieldImage;
    NiftiImage targetImage;
    
    void initImages (nifti_image *targetImage);
    
public:
    DeformationField (nifti_image *targetImage, const AffineMatrix &affine);
    DeformationField (nifti_image *targetImage, nifti_image *transformationImage);
    
    ~DeformationField ()
    {
        nifti_image_free(deformationFieldImage);
    }
    
    NiftiImage getFieldImage () const { return deformationFieldImage; }
    
    NiftiImage resampleImage (nifti_image *sourceImage, const int interpolation) const;
};

#endif
