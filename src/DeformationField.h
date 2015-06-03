#ifndef _DEFORMATION_FIELD_H_
#define _DEFORMATION_FIELD_H_

#include "config.h"
#include "nifti1_io.h"
#include "_reg_localTransformation.h"
#include "_reg_globalTransformation.h"
#include "_reg_resampling.h"
#include "nifti_image.h"
#include "aladin.h"

class DeformationField
{
protected:
    nifti_image *deformationFieldImage;
    nifti_image *targetImage;
    
    void initImages (nifti_image *targetImage)
    {
        this->targetImage = targetImage;
        
        // Create a deformation field
        nifti_image *deformationFieldImage = nifti_copy_nim_info(targetImage);
        deformationFieldImage->dim[0] = deformationFieldImage->ndim = 5;
        deformationFieldImage->dim[1] = deformationFieldImage->nx = targetImage->nx;
        deformationFieldImage->dim[2] = deformationFieldImage->ny = targetImage->ny;
        deformationFieldImage->dim[3] = deformationFieldImage->nz = targetImage->nz;
        deformationFieldImage->dim[4] = deformationFieldImage->nt = 1;
        deformationFieldImage->pixdim[4] = deformationFieldImage->dt = 1.0;
        deformationFieldImage->dim[5] = deformationFieldImage->nu = (targetImage->nz>1 ? 3 : 2);
        deformationFieldImage->dim[6] = deformationFieldImage->nv = 1;
        deformationFieldImage->dim[7] = deformationFieldImage->nw = 1;
        deformationFieldImage->nvox = size_t(deformationFieldImage->nx *
            deformationFieldImage->ny * deformationFieldImage->nz *
            deformationFieldImage->nt * deformationFieldImage->nu);
        deformationFieldImage->scl_slope = 1.0f;
        deformationFieldImage->scl_inter = 0.0f;
        
        deformationFieldImage->datatype = (sizeof(PRECISION_TYPE)==4 ? NIFTI_TYPE_FLOAT32 : NIFTI_TYPE_FLOAT64);
        deformationFieldImage->nbyper = sizeof(PRECISION_TYPE);
        deformationFieldImage->data = (void *) calloc(deformationFieldImage->nvox, deformationFieldImage->nbyper);

        // Initialise the deformation field with an identity transformation
        reg_tools_multiplyValueToImage(deformationFieldImage, deformationFieldImage, 0.0f);
        reg_getDeformationFromDisplacement(deformationFieldImage);
        deformationFieldImage->intent_p1 = DEF_FIELD;
    }
    
public:
    DeformationField ()
        : deformationFieldImage(NULL) {}
    
    DeformationField (nifti_image *targetImage, const AffineMatrix &affine)
    {
        initImages(targetImage);
        mat44 affineMatrix = affine;
        reg_affine_getDeformationField(&affineMatrix, deformationFieldImage, false, NULL);
    }
    
    DeformationField (nifti_image *targetImage, nifti_image *transformationImage)
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
                nifti_image *tempFlowField = copyCompleteImage(deformationFieldImage);
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
    
    ~DeformationField ()
    {
        nifti_image_free(deformationFieldImage);
    }
    
    nifti_image * getFieldImage () const { return deformationFieldImage; }
    
    nifti_image * resampleImage (nifti_image * sourceImage, const int interpolation) const
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
    
        return resultImage;
    }
};

#endif
