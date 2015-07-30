#ifndef _DEFORMATION_FIELD_H_
#define _DEFORMATION_FIELD_H_

#include <RcppEigen.h>

#include "config.h"
#include "nifti1_io.h"
#include "_reg_localTrans.h"
#include "_reg_globalTrans.h"
#include "_reg_resampling.h"
#include "NiftiImage.h"
#include "AffineMatrix.h"

class DeformationField
{
protected:
    NiftiImage deformationFieldImage;
    NiftiImage targetImage;
    
    void initImages (const NiftiImage &targetImage);
    
public:
    DeformationField () {}
    DeformationField (const NiftiImage &targetImage, const AffineMatrix &affine, const bool compose = false);
    DeformationField (const NiftiImage &targetImage, const NiftiImage &transformationImage, const bool compose = false);
    
    NiftiImage getFieldImage () const { return deformationFieldImage; }
    
    NiftiImage getJacobian () const;
    
    NiftiImage resampleImage (const NiftiImage &sourceImage, const int interpolation) const;
    
    template <int Dim>
    Rcpp::NumericVector findPoint (const NiftiImage &sourceImage, const Eigen::Matrix<double,Dim,1> &sourceLoc, const bool nearest) const;
    
    void compose (const DeformationField &otherField);
};

#endif
