#ifndef _DEFORMATION_FIELD_H_
#define _DEFORMATION_FIELD_H_

#include <RcppEigen.h>

#include "RNifti.h"
#include "_reg_localTrans.h"
#include "_reg_globalTrans.h"
#include "_reg_resampling.h"
#include "AffineMatrix.h"

template <typename PrecisionType>
class DeformationField
{
protected:
    RNifti::NiftiImage deformationFieldImage;
    RNifti::NiftiImage targetImage;
    std::vector<double> deformationData;
    size_t nVoxels;
    
    void initImages (const RNifti::NiftiImage &targetImage);
    void updateData ()
    {
        deformationData = deformationFieldImage.getData<double>();
        nVoxels = deformationFieldImage->nx * deformationFieldImage->ny * deformationFieldImage->nz;
    }
    
public:
    DeformationField () {}
    DeformationField (const RNifti::NiftiImage &targetImage, const AffineMatrix &affine, const bool compose = false);
    DeformationField (const RNifti::NiftiImage &targetImage, RNifti::NiftiImage &transformationImage, const bool compose = false);
    
    RNifti::NiftiImage getFieldImage () const { return deformationFieldImage; }
    
    RNifti::NiftiImage getJacobian ();
    
    RNifti::NiftiImage resampleImage (RNifti::NiftiImage &sourceImage, const int interpolation);
    
    template <int Dim>
    Rcpp::NumericVector findPoint (const RNifti::NiftiImage &sourceImage, const Eigen::Matrix<double,Dim,1> &sourceLoc, const bool nearest, const Eigen::Matrix<double,Dim,1> &start) const;
    
    void compose (const DeformationField &otherField);
};

#endif
