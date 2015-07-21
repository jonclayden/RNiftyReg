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
    
    AffineMatrix (const mat44 &matrix, const bool attributes = true);
    AffineMatrix (const nifti_image *sourceImage, const nifti_image *targetImage);
    
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
