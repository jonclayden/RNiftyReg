#ifndef _AFFINE_MATRIX_H_
#define _AFFINE_MATRIX_H_

#include "nifti1_io.h"
#include "_reg_maths.h"

class AffineMatrix : public Rcpp::NumericMatrix
{
private:
    void addAttributes ()
    {
        this->attr("class") = "affine";
        this->attr("affineType") = "niftyreg";
    }
    
public:
    AffineMatrix ()
        : Rcpp::NumericMatrix(4,4)
    {
        std::fill(this->begin(), this->end(), 0.0);
    }
    
    AffineMatrix (SEXP object)
        : Rcpp::NumericMatrix(object)
    {
        if (this->cols() != 4 || this->rows() != 4)
            throw std::runtime_error("Specified affine matrix does not have dimensions of 4x4");
    }
    
    AffineMatrix (const mat44 &matrix)
        : Rcpp::NumericMatrix(4,4)
    {
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
                (*this)(i,j) = static_cast<double>(matrix.m[i][j]);
        }
        
        addAttributes();
    }
    
    AffineMatrix (const nifti_image *sourceImage, const nifti_image *targetImage)
        : Rcpp::NumericMatrix(4,4)
    {
        std::fill(this->begin(), this->end(), 0.0);
        (*this)(0,0) = (*this)(1,1) = (*this)(2,2) = (*this)(3,3) = 1.0;
        
        mat44 *sourceMatrix, *targetMatrix;
        float sourceCentre[3], targetCentre[3], sourceRealPosition[3], targetRealPosition[3];
        
        if (sourceImage->sform_code>0)
            sourceMatrix = const_cast<mat44*>(&(sourceImage->sto_xyz));
        else
            sourceMatrix = const_cast<mat44*>(&(sourceImage->qto_xyz));
        if (targetImage->sform_code>0)
            targetMatrix = const_cast<mat44*>(&(targetImage->sto_xyz));
        else
            targetMatrix = const_cast<mat44*>(&(targetImage->qto_xyz));
        
        sourceCentre[0] = (float) (sourceImage->nx) / 2.0f;
        sourceCentre[1] = (float) (sourceImage->ny) / 2.0f;
        sourceCentre[2] = (float) (sourceImage->nz) / 2.0f;
        
        targetCentre[0] = (float) (targetImage->nx) / 2.0f;
        targetCentre[1] = (float) (targetImage->ny) / 2.0f;
        targetCentre[2] = (float) (targetImage->nz) / 2.0f;
        
        reg_mat44_mul(sourceMatrix, sourceCentre, sourceRealPosition);
        reg_mat44_mul(targetMatrix, targetCentre, targetRealPosition);
        
        // Use origins to initialise translation elements
        (*this)(0,3) = sourceRealPosition[0] - targetRealPosition[0];
        (*this)(1,3) = sourceRealPosition[1] - targetRealPosition[1];
        (*this)(2,3) = sourceRealPosition[2] - targetRealPosition[2];
        
        addAttributes();
    }
    
    operator mat44 () const
    {
        mat44 matrix;
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<4; j++)
                matrix.m[i][j] = static_cast<float>((*this)(i,j));
        }
        return matrix;
    }
    
    bool isValid () const { return ((*this)(3,3) != 0.0); }
};

#endif
