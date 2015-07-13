#include <RcppEigen.h>

#include "_reg_localTransformation.h"
#include "_reg_globalTransformation.h"
#include "_reg_resampling.h"

#include "config.h"
#include "DeformationField.h"

void DeformationField::initImages (nifti_image *targetImage)
{
    this->targetImage = NiftiImage(targetImage);
    
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

DeformationField::DeformationField (nifti_image *targetImage, const AffineMatrix &affine)
{
    initImages(targetImage);
    mat44 affineMatrix = affine;
    reg_affine_getDeformationField(&affineMatrix, deformationFieldImage, false, NULL);
}

DeformationField::DeformationField (nifti_image *targetImage, nifti_image *transformationImage)
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

NiftiImage DeformationField::resampleImage (nifti_image *sourceImage, const int interpolation) const
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
