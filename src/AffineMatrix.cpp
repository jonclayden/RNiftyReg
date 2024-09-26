#include <RcppEigen.h>

#include "AffineMatrix.h"

AffineMatrix::AffineMatrix (const mat44 &matrix, const bool attributes)
    : Rcpp::NumericMatrix(4,4)
{
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            (*this)(i,j) = static_cast<double>(matrix.m[i][j]);
    }
    
    if (attributes)
        addAttributes();
}

AffineMatrix::AffineMatrix (const Eigen::MatrixXd &matrix, const bool attributes)
    : Rcpp::NumericMatrix(4,4)
{
    if (matrix.rows() != 4 || matrix.cols() != 4)
        throw std::runtime_error("The specified matrix is not 4x4");
    
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
            (*this)(i,j) = matrix(i,j);
    }
    
    if (attributes)
        addAttributes();
}

AffineMatrix::AffineMatrix (const RNifti::NiftiImage &sourceImage, const RNifti::NiftiImage &targetImage)
    : Rcpp::NumericMatrix(4,4)
{
    std::fill(this->begin(), this->end(), 0.0);
    (*this)(0,0) = (*this)(1,1) = (*this)(2,2) = (*this)(3,3) = 1.0;
    
    const mat44 sourceMatrix = sourceImage.xform(false);
    const mat44 targetMatrix = targetImage.xform(false);
    float sourceCentre[3], targetCentre[3], sourceRealPosition[3], targetRealPosition[3];
    
    sourceCentre[0] = (float) (sourceImage->nx) / 2.0f;
    sourceCentre[1] = (float) (sourceImage->ny) / 2.0f;
    sourceCentre[2] = (float) (sourceImage->nz) / 2.0f;
    
    targetCentre[0] = (float) (targetImage->nx) / 2.0f;
    targetCentre[1] = (float) (targetImage->ny) / 2.0f;
    targetCentre[2] = (float) (targetImage->nz) / 2.0f;
    
    reg_mat44_mul(&sourceMatrix, sourceCentre, sourceRealPosition);
    reg_mat44_mul(&targetMatrix, targetCentre, targetRealPosition);
    
    // Use origins to initialise translation elements
    (*this)(0,3) = sourceRealPosition[0] - targetRealPosition[0];
    (*this)(1,3) = sourceRealPosition[1] - targetRealPosition[1];
    (*this)(2,3) = sourceRealPosition[2] - targetRealPosition[2];
    
    addAttributes();
}
