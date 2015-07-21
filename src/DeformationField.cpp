#include <RcppEigen.h>

#include "_reg_localTransformation.h"
#include "_reg_globalTransformation.h"
#include "_reg_resampling.h"

#include "config.h"
#include "DeformationField.h"

void DeformationField::initImages (const NiftiImage &targetImage)
{
    this->targetImage = targetImage;
    
    // Create a deformation field
    nifti_image *deformationField = nifti_copy_nim_info(targetImage);
    deformationField->dim[0] = deformationField->ndim = 5;
    deformationField->dim[1] = deformationField->nx = targetImage->nx;
    deformationField->dim[2] = deformationField->ny = targetImage->ny;
    deformationField->dim[3] = deformationField->nz = targetImage->nz;
    deformationField->dim[4] = deformationField->nt = 1;
    deformationField->pixdim[4] = deformationField->dt = 1.0;
    deformationField->dim[5] = deformationField->nu = (targetImage->nz>1 ? 3 : 2);
    deformationField->dim[6] = deformationField->nv = 1;
    deformationField->dim[7] = deformationField->nw = 1;
    deformationField->nvox = size_t(deformationField->nx *
        deformationField->ny * deformationField->nz *
        deformationField->nt * deformationField->nu);
    deformationField->scl_slope = 1.0f;
    deformationField->scl_inter = 0.0f;
    
    deformationField->datatype = (sizeof(PRECISION_TYPE)==4 ? NIFTI_TYPE_FLOAT32 : NIFTI_TYPE_FLOAT64);
    deformationField->nbyper = sizeof(PRECISION_TYPE);
    deformationField->data = (void *) calloc(deformationField->nvox, deformationField->nbyper);

    // Initialise the deformation field with an identity transformation
    reg_tools_multiplyValueToImage(deformationField, deformationField, 0.0f);
    reg_getDeformationFromDisplacement(deformationField);
    deformationField->intent_p1 = DEF_FIELD;
    
    this->deformationFieldImage = NiftiImage(deformationField);
}

DeformationField::DeformationField (const NiftiImage &targetImage, const AffineMatrix &affine)
{
    initImages(targetImage);
    mat44 affineMatrix = affine;
    reg_affine_getDeformationField(&affineMatrix, deformationFieldImage, false, NULL);
}

DeformationField::DeformationField (const NiftiImage &targetImage, const NiftiImage &transformationImage)
{
    initImages(targetImage);
    reg_checkAndCorrectDimension(transformationImage);
    
    switch (static_cast<int>(transformationImage->intent_p1))
    {
        case SPLINE_GRID:
        reg_spline_getDeformationField(transformationImage, deformationFieldImage, NULL, false, true);
        break;
        
        case DISP_VEL_FIELD:
        reg_getDeformationFromDisplacement(transformationImage);
        case DEF_VEL_FIELD:
        {
            nifti_image *tempFlowField = deformationFieldImage;
            reg_defField_compose(transformationImage, tempFlowField, NULL);
            tempFlowField->intent_p1 = transformationImage->intent_p1;
            tempFlowField->intent_p2 = transformationImage->intent_p2;
            reg_defField_getDeformationFieldFromFlowField(tempFlowField, deformationFieldImage, false);
            nifti_image_free(tempFlowField);
        }
        break;
        
        case SPLINE_VEL_GRID:
        reg_spline_getDefFieldFromVelocityGrid(transformationImage, deformationFieldImage, false);
        break;
        
        case DISP_FIELD:
        reg_getDeformationFromDisplacement(transformationImage);
        default:
        reg_defField_compose(transformationImage, deformationFieldImage, NULL);
        break;
    }
}

NiftiImage DeformationField::getJacobian () const
{
    // Allocate Jacobian determinant image
    nifti_image *jacobianImage = nifti_copy_nim_info(targetImage);
    jacobianImage->cal_min = 0;
    jacobianImage->cal_max = 0;
    jacobianImage->scl_slope = 1.0;
    jacobianImage->scl_inter = 0.0;
    jacobianImage->datatype = NIFTI_TYPE_FLOAT64;
    jacobianImage->nbyper = 8;
    jacobianImage->data = (void *) calloc(jacobianImage->nvox, jacobianImage->nbyper);
    
    // Calculate Jacobian determinant map
    reg_defField_getJacobianMap(deformationFieldImage, jacobianImage);
    
    return NiftiImage(jacobianImage);
}

NiftiImage DeformationField::resampleImage (const NiftiImage &sourceImage, const int interpolation) const
{
    // Allocate result image
    nifti_image *resultImage = nifti_copy_nim_info(targetImage);
    resultImage->dim[0] = resultImage->ndim = sourceImage->dim[0];
    resultImage->dim[4] = resultImage->nt = sourceImage->dim[4];
    resultImage->cal_min = sourceImage->cal_min;
    resultImage->cal_max = sourceImage->cal_max;
    resultImage->scl_slope = sourceImage->scl_slope;
    resultImage->scl_inter = sourceImage->scl_inter;
    resultImage->datatype = sourceImage->datatype;
    resultImage->nbyper = sourceImage->nbyper;
    resultImage->nvox = size_t(resultImage->dim[1]) * size_t(resultImage->dim[2]) * size_t(resultImage->dim[3]) * size_t(resultImage->dim[4]);
    resultImage->data = (void *) calloc(resultImage->nvox, resultImage->nbyper);

    // Resample source image to target space
    reg_resampleImage(sourceImage, resultImage, deformationFieldImage, NULL, interpolation, 0);

    return NiftiImage(resultImage);
}

template <int Dim>
Rcpp::NumericVector DeformationField::findPoint (const NiftiImage &sourceImage, const Eigen::Matrix<double,Dim,1> &sourceLoc, const bool nearest) const
{
    typedef Eigen::Matrix<double,Dim,1> Point;
    Point closestLoc;
    
    const PRECISION_TYPE *deformationPointer = (const PRECISION_TYPE *) deformationFieldImage->data;
    const size_t nVoxels = deformationFieldImage->nx * deformationFieldImage->ny * deformationFieldImage->nz;
    double closestDistance = R_PosInf;
    size_t closestVoxel = 0;
    
    for (size_t v=0; v<nVoxels; v++)
    {
        Point loc;
        for (int i=0; i<Dim; i++)
            loc[i] = deformationPointer[v + i*nVoxels];
        
        const double currentDistance = (loc - sourceLoc).norm();
        if (currentDistance < closestDistance)
        {
            closestDistance = currentDistance;
            closestVoxel = v;
            closestLoc = loc;
        }
    }
    
    std::vector<size_t> strides(Dim);
    strides[0] = 1;
    for (int i=1; i<Dim; i++)
        strides[i] = strides[i-1] * std::abs(deformationFieldImage->dim[i]);
    
    if (nearest || closestDistance == 0.0)
    {
        Rcpp::NumericVector result(Dim);
        result[0] = closestVoxel % deformationFieldImage->dim[1] + 1.0;
        for (int i=1; i<Dim; i++)
            result[i] = (closestVoxel / strides[i]) % deformationFieldImage->dim[i+1] + 1.0;
        
        return result;
    }
    else
    {
        Rcpp::NumericVector result(int(R_pow_di(4.0,Dim)) * 2 * Dim);
        const mat44 &xform = sourceImage.xform();
        Point offset = sourceLoc - closestLoc;
        size_t v;
        
        for (int i=0; i<Dim; i++)
            offset[i] = (offset[i] * xform.m[i][i] >= 0.0 ? 0.0 : -1.0);
        
        for (int i=0; i<4; i++)
        {
            const int xShift = i + offset[0] - 1;

            for (int j=0; j<4; j++)
            {
                const int yShift = j + offset[1] - 1;
                
                if (Dim == 2)
                {
                    const size_t v = closestVoxel + xShift + yShift * strides[1];
                    const size_t w = size_t(i + j*4);
                
                    result[2*Dim*w] = deformationPointer[v];
                    result[2*Dim*w + 1] = deformationPointer[v + nVoxels];
                    result[2*Dim*w + 2] = v % deformationFieldImage->dim[1] + 1.0;
                    result[2*Dim*w + 3] = (v / strides[1]) % deformationFieldImage->dim[2] + 1.0;
                }
                else
                {
                    for (int k=0; k<4; k++)
                    {
                        const int zShift = k + offset[2] - 1;
                        const size_t v = closestVoxel + xShift + yShift * strides[1] + zShift * strides[2];
                        const size_t w = size_t(i + j*4 + k*16);
                        
                        result[2*Dim*w] = deformationPointer[v];
                        result[2*Dim*w + 1] = deformationPointer[v + nVoxels];
                        result[2*Dim*w + 2] = deformationPointer[v + 2*nVoxels];
                        result[2*Dim*w + 3] = v % deformationFieldImage->dim[1] + 1.0;
                        result[2*Dim*w + 4] = (v / strides[1]) % deformationFieldImage->dim[2] + 1.0;
                        result[2*Dim*w + 5] = (v / strides[2]) % deformationFieldImage->dim[3] + 1.0;
                    }
                }
            }
        }
        
        return result;
    }
}

template
Rcpp::NumericVector DeformationField::findPoint (const NiftiImage &sourceImage, const Eigen::Matrix<double,2,1> &sourceLoc, const bool nearest) const;

template
Rcpp::NumericVector DeformationField::findPoint (const NiftiImage &sourceImage, const Eigen::Matrix<double,3,1> &sourceLoc, const bool nearest) const;
