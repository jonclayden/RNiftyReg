// Registration types (degrees of freedom) - these constants have to have these values
#define RIGID 0
#define AFFINE 1

// Precision levels for final interpolation
#define INTERP_PREC_SOURCE 0
#define INTERP_PREC_DOUBLE 1

// Working precision for the source and target images: float and double are valid for the NiftyReg code
#define PRECISION_TYPE double

// Convergence criterion
#define CONVERGENCE_EPS 0.00001

#define JH_TRI 0
#define JH_PARZEN_WIN 1
#define JH_PW_APPROX 2

#define FOLDING_CORRECTION_STEP 20

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_blockMatching.h"
#include "_reg_tools.h"
#include "_reg_bspline.h"
#include "_reg_bspline_comp.h"
#include "_reg_mutualinformation.h"
#include "_reg_ssd.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "niftyreg.h"

extern "C"
SEXP reg_aladin (SEXP source, SEXP target, SEXP type, SEXP finalPrecision, SEXP nLevels, SEXP maxIterations, SEXP useBlockPercentage, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP verbose)
{
    int i, j, levels = *(INTEGER(nLevels));
    SEXP returnValue, data, completedIterations;
    
    bool affineProvided = !isNull(affineComponents);
    
    nifti_image *sourceImage = s4_image_to_struct(source);
    nifti_image *targetImage = s4_image_to_struct(target);
    nifti_image *targetMaskImage = NULL;
    mat44 *affineTransformation = NULL;
    
    if (!isNull(targetMask))
        targetMaskImage = s4_image_to_struct(targetMask);
    
    if (affineProvided)
    {
        affineTransformation = (mat44 *) calloc(1, sizeof(mat44));
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
                affineTransformation->m[i][j] = (float) REAL(affineComponents)[(j*4)+i];
        }
    }
    
    int regType = (strcmp(CHAR(STRING_ELT(type,0)),"rigid")==0 ? RIGID : AFFINE);
    int precisionType = (strcmp(CHAR(STRING_ELT(finalPrecision,0)),"source")==0 ? INTERP_PREC_SOURCE : INTERP_PREC_DOUBLE);
    
    aladin_result result = do_reg_aladin(sourceImage, targetImage, regType, precisionType, *(INTEGER(nLevels)), *(INTEGER(maxIterations)), *(INTEGER(useBlockPercentage)), *(INTEGER(finalInterpolation)), targetMaskImage, affineTransformation, (*(INTEGER(verbose)) == 1));
    
    PROTECT(returnValue = NEW_LIST(3));
    
    if (result.image->datatype == DT_INT32)
    {
        // Integer-valued data went in, and precision must be "source"
        PROTECT(data = NEW_INTEGER((R_len_t) result.image->nvox));
        for (size_t i = 0; i < result.image->nvox; i++)
            INTEGER(data)[i] = ((int *) result.image->data)[i];
    }
    else
    {
        // All other cases
        PROTECT(data = NEW_NUMERIC((R_len_t) result.image->nvox));
        for (size_t i = 0; i < result.image->nvox; i++)
            REAL(data)[i] = ((double *) result.image->data)[i];
    }
    
    SET_ELEMENT(returnValue, 0, data);
    
    if (!affineProvided)
        PROTECT(affineComponents = NEW_NUMERIC(16));
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
            REAL(affineComponents)[(j*4)+i] = (double) result.affine->m[i][j];
    }
    
    SET_ELEMENT(returnValue, 1, affineComponents);
    
    PROTECT(completedIterations = NEW_INTEGER(levels));
    for (i = 0; i < levels; i++)
        INTEGER(completedIterations)[i] = result.completedIterations[i];
    
    SET_ELEMENT(returnValue, 2, completedIterations);
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    if (targetMaskImage != NULL)
        nifti_image_free(targetMaskImage);
    nifti_image_free(result.image);
    free(result.affine);
    
    UNPROTECT(affineProvided ? 3 : 4);
    
    return returnValue;
}

extern "C"
SEXP reg_f3d (SEXP source, SEXP target, SEXP finalPrecision, SEXP nLevels, SEXP maxIterations, SEXP nBins, SEXP bendingEnergyWeight, SEXP jacobianWeight, SEXP finalSpacing, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP initControl, SEXP verbose)
{
    int i, j, levels = *(INTEGER(nLevels));
    SEXP returnValue, data, controlPoints, completedIterations;
    
    bool affineProvided = !isNull(affineComponents);
    
    nifti_image *sourceImage = s4_image_to_struct(source);
    nifti_image *targetImage = s4_image_to_struct(target);
    nifti_image *targetMaskImage = NULL;
    nifti_image *controlPointImage = NULL;
    mat44 *affineTransformation = NULL;
    
    if (!isNull(targetMask))
        targetMaskImage = s4_image_to_struct(targetMask);
    if (!isNull(initControl))
        controlPointImage = s4_image_to_struct(initControl);
    
    if (affineProvided)
    {
        affineTransformation = (mat44 *) calloc(1, sizeof(mat44));
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
                affineTransformation->m[i][j] = (float) REAL(affineComponents)[(j*4)+i];
        }
    }
    else if (isNull(initControl))
    {
        affineTransformation = (mat44 *) calloc(1, sizeof(mat44));
        for (i = 0; i < 4; i++)
        {
            for (j = 0; j < 4; j++)
                affineTransformation->m[i][j] = (i == j ? 1.0 : 0.0);
        }
    }
    
    float spacing[3];
    for (i = 0; i < 3; i++)
        spacing[i] = (float) REAL(finalSpacing)[i];
    
    int precisionType = (strcmp(CHAR(STRING_ELT(finalPrecision,0)),"source")==0 ? INTERP_PREC_SOURCE : INTERP_PREC_DOUBLE);
    
    f3d_result result = do_reg_f3d(sourceImage, targetImage, precisionType, *(INTEGER(nLevels)), *(INTEGER(maxIterations)), *(INTEGER(finalInterpolation)), targetMaskImage, controlPointImage, affineTransformation, *(INTEGER(nBins)), spacing, (float) *(REAL(bendingEnergyWeight)), (float) *(REAL(jacobianWeight)), (*(INTEGER(verbose)) == 1));
    
    PROTECT(returnValue = NEW_LIST(3));
    
    if (result.image->datatype == DT_INT32)
    {
        // Integer-valued data went in, and precision must be "source"
        PROTECT(data = NEW_INTEGER((R_len_t) result.image->nvox));
        for (size_t i = 0; i < result.image->nvox; i++)
            INTEGER(data)[i] = ((int *) result.image->data)[i];
    }
    else
    {
        // All other cases
        PROTECT(data = NEW_NUMERIC((R_len_t) result.image->nvox));
        for (size_t i = 0; i < result.image->nvox; i++)
            REAL(data)[i] = ((double *) result.image->data)[i];
    }
    
    SET_ELEMENT(returnValue, 0, data);
    
    PROTECT(controlPoints = NEW_NUMERIC((R_len_t) result.controlPoints->nvox));
    for (size_t i = 0; i < result.controlPoints->nvox; i++)
        REAL(controlPoints)[i] = ((double *) result.controlPoints->data)[i];
    
    SET_ELEMENT(returnValue, 1, controlPoints);
    
    PROTECT(completedIterations = NEW_INTEGER(levels));
    for (i = 0; i < levels; i++)
        INTEGER(completedIterations)[i] = result.completedIterations[i];
    
    SET_ELEMENT(returnValue, 2, completedIterations);
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    if (targetMaskImage != NULL)
        nifti_image_free(targetMaskImage);
    nifti_image_free(result.image);
    nifti_image_free(result.controlPoints);
    if (affineProvided || isNull(initControl))
        free(affineTransformation);
    
    UNPROTECT(4);
    
    return returnValue;
}

// Convert an S4 "nifti" object, as defined in the oro.nifti package, to a "nifti_image" struct
nifti_image * s4_image_to_struct (SEXP object)
{
    int i;
    nifti_1_header header;
    
    header.sizeof_hdr = 348;
    
    for (i = 0; i < 8; i++)
        header.dim[i] = (short) INTEGER(GET_SLOT(object, install("dim_")))[i];
    
    header.intent_p1 = (float) *(REAL(GET_SLOT(object, install("intent_p1"))));
    header.intent_p2 = (float) *(REAL(GET_SLOT(object, install("intent_p2"))));
    header.intent_p3 = (float) *(REAL(GET_SLOT(object, install("intent_p3"))));
    header.intent_code = (short) *(INTEGER(GET_SLOT(object, install("intent_code"))));
    
    header.datatype = (short) *(INTEGER(GET_SLOT(object, install("datatype"))));
    header.bitpix = (short) *(INTEGER(GET_SLOT(object, install("bitpix"))));
    
    header.slice_start = (short) *(INTEGER(GET_SLOT(object, install("slice_start"))));
    header.slice_end = (short) *(INTEGER(GET_SLOT(object, install("slice_end"))));
    header.slice_code = (char) *(INTEGER(GET_SLOT(object, install("slice_code"))));
    header.slice_duration = (float) *(REAL(GET_SLOT(object, install("slice_duration"))));
    
    for (i = 0; i < 8; i++) 
        header.pixdim[i] = (float) REAL(GET_SLOT(object, install("pixdim")))[i];
    header.xyzt_units = (char) *(INTEGER(GET_SLOT(object, install("xyzt_units"))));
    
    header.vox_offset = (float) *(REAL(GET_SLOT(object, install("vox_offset"))));
    
    header.scl_slope = (float) *(REAL(GET_SLOT(object, install("scl_slope"))));
    header.scl_inter = (float) *(REAL(GET_SLOT(object, install("scl_inter"))));
    header.toffset = (float) *(REAL(GET_SLOT(object, install("toffset"))));
    
    header.cal_max = (float) *(REAL(GET_SLOT(object, install("cal_max"))));
    header.cal_min = (float) *(REAL(GET_SLOT(object, install("cal_min"))));
    header.glmax = header.glmin = 0;
    
    strcpy(header.descrip, CHAR(STRING_ELT(GET_SLOT(object, install("descrip")),0)));
    strcpy(header.aux_file, CHAR(STRING_ELT(GET_SLOT(object, install("aux_file")),0)));
    strcpy(header.intent_name, CHAR(STRING_ELT(GET_SLOT(object, install("intent_name")),0)));
    strcpy(header.magic, CHAR(STRING_ELT(GET_SLOT(object, install("magic")),0)));
    
    header.qform_code = (short) *(INTEGER(GET_SLOT(object, install("qform_code"))));
    header.sform_code = (short) *(INTEGER(GET_SLOT(object, install("sform_code"))));
    
    header.quatern_b = (float) *(REAL(GET_SLOT(object, install("quatern_b"))));
    header.quatern_c = (float) *(REAL(GET_SLOT(object, install("quatern_c"))));
    header.quatern_d = (float) *(REAL(GET_SLOT(object, install("quatern_d"))));
    header.qoffset_x = (float) *(REAL(GET_SLOT(object, install("qoffset_x"))));
    header.qoffset_y = (float) *(REAL(GET_SLOT(object, install("qoffset_y"))));
    header.qoffset_z = (float) *(REAL(GET_SLOT(object, install("qoffset_z"))));
    
    for (i = 0; i < 4; i++)
    {
        header.srow_x[i] = (float) REAL(GET_SLOT(object, install("srow_x")))[i];
        header.srow_y[i] = (float) REAL(GET_SLOT(object, install("srow_y")))[i];
        header.srow_z[i] = (float) REAL(GET_SLOT(object, install("srow_z")))[i];
    }
    
    if (header.datatype == DT_UINT8 || header.datatype == DT_INT16 || header.datatype == DT_INT32 || header.datatype == DT_INT8 || header.datatype == DT_UINT16 || header.datatype == DT_UINT32)
        header.datatype = DT_INT32;
    else if (header.datatype == DT_FLOAT32 || header.datatype == DT_FLOAT64)
        header.datatype = DT_FLOAT64;  // This assumes that sizeof(double) == 8
    else
    {
        error("Data type %d is not supported", header.datatype);
        return NULL;
    }
    
    nifti_image *image = nifti_convert_nhdr2nim(header, NULL);
    
    size_t dataSize = nifti_get_volsize(image);
    image->data = calloc(1, dataSize);
    if (header.datatype == DT_INT32)
        memcpy(image->data, INTEGER(GET_SLOT(object, install(".Data"))), dataSize);
    else
        memcpy(image->data, REAL(GET_SLOT(object, install(".Data"))), dataSize);
    
    return image;
}

// Does the specified update matrix represent convergence?
bool reg_test_convergence (mat44 *updateMatrix)
{
    bool convergence = true;
    float referenceValue;
    
    // The fourth row is always [0 0 0 1] so we don't test it
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // The diagonal represents scale factors, so 1 indicates no change
            referenceValue = (i==j ? 1.0f : 0.0f);
            if ((fabsf(updateMatrix->m[i][j]) - referenceValue) > CONVERGENCE_EPS)
                convergence = false;
        }
    }

    return convergence;
}

nifti_image * copy_complete_nifti_image (nifti_image *source)
{
    size_t dataSize = nifti_get_volsize(source);
    nifti_image *destination = nifti_copy_nim_info(source);
    destination->data = calloc(1, dataSize);
    memcpy(destination->data, source->data, dataSize);
    
    return destination;
}

nifti_image * create_position_field (nifti_image *templateImage, bool twoDimRegistration)
{
    nifti_image *positionFieldImage = nifti_copy_nim_info(templateImage);
    
    positionFieldImage->dim[0] = positionFieldImage->ndim = 5;
    positionFieldImage->dim[1] = positionFieldImage->nx = templateImage->nx;
    positionFieldImage->dim[2] = positionFieldImage->ny = templateImage->ny;
    positionFieldImage->dim[3] = positionFieldImage->nz = templateImage->nz;
    positionFieldImage->dim[4] = positionFieldImage->nt = 1;
    positionFieldImage->pixdim[4] = positionFieldImage->dt = 1.0;
    positionFieldImage->dim[5] = positionFieldImage->nu = (twoDimRegistration ? 2 : 3);
    positionFieldImage->pixdim[5] = positionFieldImage->du = 1.0;
    positionFieldImage->dim[6] = positionFieldImage->nv = 1;
    positionFieldImage->pixdim[6] = positionFieldImage->dv = 1.0;
    positionFieldImage->dim[7] = positionFieldImage->nw = 1;
    positionFieldImage->pixdim[7] = positionFieldImage->dw = 1.0;
    positionFieldImage->nvox = positionFieldImage->nx * positionFieldImage->ny * positionFieldImage->nz * positionFieldImage->nt * positionFieldImage->nu;
    positionFieldImage->datatype = (sizeof(PRECISION_TYPE)==4 ? NIFTI_TYPE_FLOAT32 : NIFTI_TYPE_FLOAT64);
    positionFieldImage->nbyper = sizeof(PRECISION_TYPE);
    positionFieldImage->data = calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
    
    return positionFieldImage;
}

// Run the "aladin" registration algorithm
aladin_result do_reg_aladin (nifti_image *sourceImage, nifti_image *targetImage, int type, int finalPrecision, int nLevels, int maxIterations, int useBlockPercentage, int finalInterpolation, nifti_image *targetMaskImage, mat44 *affineTransformation, bool verbose)
{
    int i;
    bool usingTargetMask = (targetMaskImage != NULL);
    bool twoDimRegistration = (sourceImage->nz == 1 || targetImage->nz == 1);
    float sourceBGValue = 0.0;
    nifti_image *resultImage, *positionFieldImage = NULL;
    
    int *completedIterations = (int *) calloc(nLevels, sizeof(int));
    
    // Initial affine matrix is the identity
    if (affineTransformation == NULL)
    {
        affineTransformation = (mat44 *) calloc(1, sizeof(mat44));
        for (i = 0; i < 4; i++)
            affineTransformation->m[i][i] = 1.0;

        mat44 *sourceMatrix, *targetMatrix;
        float sourceCentre[3], targetCentre[3], sourceRealPosition[3], targetRealPosition[3];

    	if (sourceImage->sform_code>0)
    		sourceMatrix = &(sourceImage->sto_xyz);
    	else
    	    sourceMatrix = &(sourceImage->qto_xyz);
    	if (targetImage->sform_code>0)
    		targetMatrix = &(targetImage->sto_xyz);
    	else
    	    targetMatrix = &(targetImage->qto_xyz);

    	sourceCentre[0] = (float) (sourceImage->nx) / 2.0f;
    	sourceCentre[1] = (float) (sourceImage->ny) / 2.0f;
    	sourceCentre[2] = (float) (sourceImage->nz) / 2.0f;

    	targetCentre[0] = (float) (targetImage->nx) / 2.0f;
    	targetCentre[1] = (float) (targetImage->ny) / 2.0f;
    	targetCentre[2] = (float) (targetImage->nz) / 2.0f;

    	reg_mat44_mul(sourceMatrix, sourceCentre, sourceRealPosition);
    	reg_mat44_mul(targetMatrix, targetCentre, targetRealPosition);

    	// Use origins to initialise translation elements
    	affineTransformation->m[0][3] = sourceRealPosition[0] - targetRealPosition[0];
    	affineTransformation->m[1][3] = sourceRealPosition[1] - targetRealPosition[1];
    	affineTransformation->m[2][3] = sourceRealPosition[2] - targetRealPosition[2];
    }
    
    // Binarise the mask image
    if (usingTargetMask)
        reg_tool_binarise_image(targetMaskImage);
    
    for (int level = 0; level < nLevels; level++)
    {
        nifti_image *sourceImageCopy = copy_complete_nifti_image(sourceImage);
        nifti_image *targetImageCopy = copy_complete_nifti_image(targetImage);
        nifti_image *targetMaskImageCopy = NULL;
        if (usingTargetMask)
            targetMaskImageCopy = copy_complete_nifti_image(targetMaskImage);
        
        reg_changeDatatype<PRECISION_TYPE>(sourceImageCopy);
        reg_changeDatatype<PRECISION_TYPE>(targetImageCopy);
        
        for (int l = level; l < nLevels-1; l++)
        {
            int ratio = (int) powf(2.0f, l+1.0f);

            bool sourceDownsampleAxis[8] = {true,true,true,true,true,true,true,true};
            if ((sourceImage->nx/ratio) < 32) sourceDownsampleAxis[1] = false;
            if ((sourceImage->ny/ratio) < 32) sourceDownsampleAxis[2] = false;
            if ((sourceImage->nz/ratio) < 32) sourceDownsampleAxis[3] = false;
            reg_downsampleImage<PRECISION_TYPE>(sourceImageCopy, 1, sourceDownsampleAxis);

            bool targetDownsampleAxis[8] = {true,true,true,true,true,true,true,true};
            if ((targetImage->nx/ratio) < 32) targetDownsampleAxis[1] = false;
            if ((targetImage->ny/ratio) < 32) targetDownsampleAxis[2] = false;
            if ((targetImage->nz/ratio) < 32) targetDownsampleAxis[3] = false;
            reg_downsampleImage<PRECISION_TYPE>(targetImageCopy, 1, targetDownsampleAxis);

            if (usingTargetMask)
                reg_downsampleImage<PRECISION_TYPE>(targetMaskImageCopy, 0, targetDownsampleAxis);
        }
        
        int activeVoxelNumber = 0;
        int *targetMask = (int *) calloc(targetImageCopy->nvox, sizeof(int));
        if (usingTargetMask)
        {
            reg_tool_binaryImage2int(targetMaskImageCopy, targetMask, activeVoxelNumber);
            nifti_image_free(targetMaskImageCopy);
        }
        else
        {
            for (unsigned int j = 0; j < targetImageCopy->nvox; j++)
                targetMask[j] = j;
            activeVoxelNumber = targetImageCopy->nvox;
        }
        
        // Allocate the deformation field image
        positionFieldImage = create_position_field(targetImageCopy, twoDimRegistration);
        
        // Allocate the result image
        resultImage = nifti_copy_nim_info(targetImageCopy);
        resultImage->datatype = sourceImageCopy->datatype;
        resultImage->nbyper = sourceImageCopy->nbyper;
        resultImage->data = calloc(1, nifti_get_volsize(resultImage));

        // Initialise the block matching - all the blocks are used during the first level
        _reg_blockMatchingParam blockMatchingParams;
        initialise_block_matching_method(targetImageCopy, &blockMatchingParams, (level==0 ? 100 : useBlockPercentage), 50, targetMask, 0);

        if (verbose)
        {
            // Display some parameters specific to the current level
            Rprintf("Current level %i / %i\n", level+1, nLevels);
            Rprintf("Target image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n", targetImageCopy->nx, targetImageCopy->ny, targetImageCopy->nz, targetImageCopy->dx, targetImageCopy->dy, targetImageCopy->dz);
            Rprintf("Source image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n", sourceImageCopy->nx, sourceImageCopy->ny, sourceImageCopy->nz, sourceImageCopy->dx, sourceImageCopy->dy, sourceImageCopy->dz);
            if (twoDimRegistration)
                Rprintf("Block size = [4 4 1]\n");
            else
                Rprintf("Block size = [4 4 4]\n");
            Rprintf("Block number = [%i %i %i]\n", blockMatchingParams.blockNumber[0], blockMatchingParams.blockNumber[1], blockMatchingParams.blockNumber[2]);
            Rprintf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
            reg_mat44_disp(affineTransformation, (char *) "Initial affine transformation");
        }

        int nLoops = ((type==AFFINE && level==0) ? 2 : 1);
        int currentType, iteration;
        mat44 updateAffineMatrix;
        
        for (i = 0; i < nLoops; i++)
        {
            // The first optimisation is rigid even if the final scope is affine
            currentType = ((i==0 && nLoops==2) ? RIGID : type);
            iteration = 0;
            
            // Twice as many iterations are performed during the first level
            while (iteration < (level==0 ? 2*maxIterations : maxIterations))
            {
                // Compute the affine transformation deformation field
                reg_affine_positionField(affineTransformation, targetImageCopy, positionFieldImage);
                
                // Resample the source image
                reg_resampleSourceImage<PRECISION_TYPE>(targetImageCopy, sourceImageCopy, resultImage, positionFieldImage, targetMask, 1, sourceBGValue);
                
                // Compute the correspondances between blocks - this is the expensive bit
                block_matching_method<PRECISION_TYPE>(targetImageCopy, resultImage, &blockMatchingParams, targetMask);
                
                // Optimise the update matrix
                optimize(&blockMatchingParams, &updateAffineMatrix, currentType);
                
                // Update the affine transformation matrix
                *affineTransformation = reg_mat44_mul(affineTransformation, &(updateAffineMatrix));
                    
                if (reg_test_convergence(&updateAffineMatrix))
                    break;
                
                iteration++;
            }
        }
        
        completedIterations[level] = iteration;

        free(targetMask);
        nifti_image_free(resultImage);
        nifti_image_free(targetImageCopy);
        nifti_image_free(sourceImageCopy);
        
        if (level < (nLevels - 1))
            nifti_image_free(positionFieldImage);
        
        if (verbose)
        {
            reg_mat44_disp(affineTransformation, (char *)"Final affine transformation");
            Rprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
        }
    }
    
    if (nLevels == 0)
        positionFieldImage = create_position_field(targetImage, twoDimRegistration);
    
    // The corresponding deformation field is evaluated and saved
    reg_affine_positionField(affineTransformation, targetImage, positionFieldImage);
    
    // The source data type is changed for precision if requested
    if (finalPrecision == INTERP_PREC_DOUBLE)
        reg_changeDatatype<double>(sourceImage);

    // The result image is resampled using a cubic spline interpolation
    resultImage = nifti_copy_nim_info(targetImage);
    resultImage->cal_min = sourceImage->cal_min;
    resultImage->cal_max = sourceImage->cal_max;
    resultImage->scl_slope = sourceImage->scl_slope;
    resultImage->scl_inter = sourceImage->scl_inter;
    resultImage->datatype = sourceImage->datatype;
    resultImage->nbyper = sourceImage->nbyper;
    resultImage->data = calloc(resultImage->nvox, resultImage->nbyper);
    reg_resampleSourceImage<PRECISION_TYPE>(targetImage, sourceImage, resultImage, positionFieldImage, NULL, finalInterpolation, sourceBGValue);
    
    nifti_image_free(positionFieldImage);
    
    aladin_result result;
    result.image = resultImage;
    result.affine = affineTransformation;
    result.completedIterations = completedIterations;
    
    return result;
}

f3d_result do_reg_f3d (nifti_image *sourceImage, nifti_image *targetImage, int finalPrecision, int nLevels, int maxIterations, int finalInterpolation, nifti_image *targetMaskImage, nifti_image *controlPointImage, mat44 *affineTransformation, int nBins, float *spacing, float bendingEnergyWeight, float jacobianWeight, bool verbose)
{
    bool usingTargetMask = (targetMaskImage != NULL);
    bool controlPointImageProvided = (controlPointImage != NULL);
    bool twoDimRegistration = (sourceImage->nz == 1 || targetImage->nz == 1);
    PRECISION_TYPE sourcePaddingValue = std::numeric_limits<float>::quiet_NaN();
    nifti_image *resultImage, *positionFieldImage = NULL;
    
    int *completedIterations = (int *) calloc(nLevels, sizeof(int));
    
    // This is due to the extrapolation of the joint histogram using the Parzen window
    nBins += 4;
    
    if (spacing[0] < 0)
        spacing[0] *= -1.0f * targetImage->dx;
    if (spacing[1] < 1)
        spacing[1] *= -1.0f * targetImage->dy;
    if (spacing[2] < 2)
        spacing[2] *= -1.0f * targetImage->dz;
    
    for (int level = 0; level < nLevels; level++)
    {
        nifti_image *sourceImageCopy = copy_complete_nifti_image(sourceImage);
        nifti_image *targetImageCopy = copy_complete_nifti_image(targetImage);
        nifti_image *targetMaskImageCopy = NULL;
        if (usingTargetMask)
            targetMaskImageCopy = copy_complete_nifti_image(targetMaskImage);
        
        reg_changeDatatype<PRECISION_TYPE>(sourceImageCopy);
        reg_changeDatatype<PRECISION_TYPE>(targetImageCopy);
        
        for (int l = level; l < nLevels-1; l++)
        {
            int ratio = (int) powf(2.0f, l+1.0f);

            bool sourceDownsampleAxis[8] = {true,true,true,true,true,true,true,true};
            if ((sourceImage->nx/ratio) < 32) sourceDownsampleAxis[1] = false;
            if ((sourceImage->ny/ratio) < 32) sourceDownsampleAxis[2] = false;
            if ((sourceImage->nz/ratio) < 32) sourceDownsampleAxis[3] = false;
            reg_downsampleImage<PRECISION_TYPE>(sourceImageCopy, 1, sourceDownsampleAxis);

            bool targetDownsampleAxis[8] = {true,true,true,true,true,true,true,true};
            if ((targetImage->nx/ratio) < 32) targetDownsampleAxis[1] = false;
            if ((targetImage->ny/ratio) < 32) targetDownsampleAxis[2] = false;
            if ((targetImage->nz/ratio) < 32) targetDownsampleAxis[3] = false;
            reg_downsampleImage<PRECISION_TYPE>(targetImageCopy, 1, targetDownsampleAxis);

            if (usingTargetMask)
                reg_downsampleImage<PRECISION_TYPE>(targetMaskImageCopy, 0, targetDownsampleAxis);
        }
        
        int activeVoxelNumber = 0;
        int *targetMask = (int *) calloc(targetImageCopy->nvox, sizeof(int));
        if (usingTargetMask)
        {
            reg_tool_binaryImage2int(targetMaskImageCopy, targetMask, activeVoxelNumber);
            nifti_image_free(targetMaskImageCopy);
        }
        else
        {
            for (unsigned int j = 0; j < targetImageCopy->nvox; j++)
                targetMask[j] = j;
            activeVoxelNumber = targetImageCopy->nvox;
        }

        // Resample the target and source images to [0,nBins-1]
        // The images are then shifted by two, which is the suport of the spline used by the parzen window filling of the joint histogram
        reg_intensityRescale(targetImageCopy, 2.0f, (float)nBins-3.0f, -FLT_MAX, FLT_MAX);
        reg_intensityRescale(sourceImageCopy, 2.0f, (float)nBins-3.0f, -FLT_MAX, FLT_MAX);

        if (level == 0)
        {
            if (!controlPointImageProvided)
            {
                // Allocate the control point image
                int dim_cpp[8];
                float gridSpacing[3];
                dim_cpp[0] = 5;
                gridSpacing[0] = spacing[0] * powf(2.0f, (float)(nLevels-1));
                dim_cpp[1] = (int) floor(targetImageCopy->nx*targetImageCopy->dx/gridSpacing[0]) + 5;
                gridSpacing[1] = spacing[1] * powf(2.0f, (float)(nLevels-1));
                dim_cpp[2] = (int) floor(targetImageCopy->ny*targetImageCopy->dy/gridSpacing[1]) + 5;
                if (twoDimRegistration)
                {
                    gridSpacing[2] = 1.0f;
                	dim_cpp[3] = 1;
                	dim_cpp[5] = 2;
                }
                else
                {
                    gridSpacing[2] = spacing[2] * powf(2.0f, (float)(nLevels-1));
                    dim_cpp[3] = (int) floor(targetImageCopy->nz*targetImageCopy->dz/gridSpacing[2]) + 5;
                	dim_cpp[5] = 3;
                }
                dim_cpp[4] = dim_cpp[6] = dim_cpp[7] = 1;
                
                if (sizeof(PRECISION_TYPE) == 4)
                    controlPointImage = nifti_make_new_nim(dim_cpp, NIFTI_TYPE_FLOAT32, true);
                else
                    controlPointImage = nifti_make_new_nim(dim_cpp, NIFTI_TYPE_FLOAT64, true);
                
                controlPointImage->cal_min = 0;
                controlPointImage->cal_max = 0;
                controlPointImage->pixdim[0] = 1.0f;
                controlPointImage->pixdim[1] = controlPointImage->dx = gridSpacing[0];
                controlPointImage->pixdim[2] = controlPointImage->dy = gridSpacing[1];
                if (twoDimRegistration)
                    controlPointImage->pixdim[3] = controlPointImage->dz = 1.0f;
                else
                    controlPointImage->pixdim[3] = controlPointImage->dz = gridSpacing[2];
                controlPointImage->pixdim[4] = controlPointImage->dt = 1.0f;
                controlPointImage->pixdim[5] = controlPointImage->du = 1.0f;
                controlPointImage->pixdim[6] = controlPointImage->dv = 1.0f;
                controlPointImage->pixdim[7] = controlPointImage->dw = 1.0f;
                controlPointImage->qform_code = targetImageCopy->qform_code;
                controlPointImage->sform_code = targetImageCopy->sform_code;
            }
        }
        else
            reg_bspline_refineControlPointGrid(targetImageCopy, controlPointImage);

		// The qform (and sform) are set for the control point position image
        float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac;
        nifti_mat44_to_quatern(targetImageCopy->qto_xyz, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac);
        controlPointImage->quatern_b = qb;
        controlPointImage->quatern_c = qc;
        controlPointImage->quatern_d = qd;
        controlPointImage->qfac = qfac;

        controlPointImage->qto_xyz = nifti_quatern_to_mat44(qb, qc, qd, qx, qy, qz, controlPointImage->dx, controlPointImage->dy, controlPointImage->dz, qfac);

        // Origin is shifted from 1 control point in the qform
        float originIndex[3];
        float originReal[3];
        originIndex[0] = -1.0f;
        originIndex[1] = -1.0f;
        originIndex[2] = 0.0f;
        if (targetImageCopy->nz > 1)
            originIndex[2] = -1.0f;
        reg_mat44_mul(&(controlPointImage->qto_xyz), originIndex, originReal);
        if (controlPointImage->qform_code == 0)
            controlPointImage->qform_code = 1;
        controlPointImage->qto_xyz.m[0][3] = controlPointImage->qoffset_x = originReal[0];
        controlPointImage->qto_xyz.m[1][3] = controlPointImage->qoffset_y = originReal[1];
        controlPointImage->qto_xyz.m[2][3] = controlPointImage->qoffset_z = originReal[2];

        controlPointImage->qto_ijk = nifti_mat44_inverse(controlPointImage->qto_xyz);

        if (controlPointImage->sform_code > 0)
        {
			nifti_mat44_to_quatern(targetImageCopy->sto_xyz, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac);

			controlPointImage->sto_xyz = nifti_quatern_to_mat44(qb, qc, qd, qx, qy, qz, controlPointImage->dx, controlPointImage->dy, controlPointImage->dz, qfac);

			// Origin is shifted from 1 control point in the sform
			originIndex[0] = -1.0f;
			originIndex[1] = -1.0f;
			originIndex[2] = 0.0f;
			if (targetImageCopy->nz > 1)
			    originIndex[2] = -1.0f;
			reg_mat44_mul(&(controlPointImage->sto_xyz), originIndex, originReal);
			controlPointImage->sto_xyz.m[0][3] = originReal[0];
			controlPointImage->sto_xyz.m[1][3] = originReal[1];
			controlPointImage->sto_xyz.m[2][3] = originReal[2];

            controlPointImage->sto_ijk = nifti_mat44_inverse(controlPointImage->sto_xyz);
        }

        if (level==0 && !controlPointImageProvided)
        {
            // The control point position image is initialised with the affine transformation
            if (reg_bspline_initialiseControlPointGridWithAffine(affineTransformation, controlPointImage))
                error("Failed to initialise control point image");
        }

        mat44 *cppMatrix_xyz;
        mat44 *targetMatrix_ijk;
        mat44 *sourceMatrix_xyz;
        if (controlPointImage->sform_code)
            cppMatrix_xyz = &(controlPointImage->sto_xyz);
        else
            cppMatrix_xyz = &(controlPointImage->qto_xyz);
        if (targetImageCopy->sform_code)
            targetMatrix_ijk = &(targetImageCopy->sto_ijk);
        else
            targetMatrix_ijk = &(targetImageCopy->qto_ijk);
        if (sourceImageCopy->sform_code)
            sourceMatrix_xyz = &(sourceImageCopy->sto_xyz);
        else
            sourceMatrix_xyz = &(sourceImageCopy->qto_xyz);

        mat33 reorient_NMI_gradient;
        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int  j=0; j<3; j++)
                reorient_NMI_gradient.m[i][j] = sourceMatrix_xyz->m[i][j];
        }

        /* allocate the deformation Field image */
        positionFieldImage = create_position_field(targetImageCopy, twoDimRegistration);

        /* allocate the result image */
        resultImage = nifti_copy_nim_info(targetImageCopy);
        resultImage->datatype = sourceImageCopy->datatype;
        resultImage->nbyper = sourceImageCopy->nbyper;
        resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);

        if (verbose)
        {
    		printf("Current level %i / %i\n", level+1, nLevels);
    		printf("Target image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n", targetImageCopy->nx, targetImageCopy->ny, targetImageCopy->nz, targetImageCopy->dx, targetImageCopy->dy, targetImageCopy->dz);
    		printf("Source image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n", sourceImageCopy->nx, sourceImageCopy->ny, sourceImageCopy->nz, sourceImageCopy->dx, sourceImageCopy->dy, sourceImageCopy->dz);
    		printf("\t%ix%ix%i control points (%i DoF)\n", controlPointImage->nx, controlPointImage->ny, controlPointImage->nz, (int) controlPointImage->nvox);
    		printf("\t%gx%gx%g mm\n", controlPointImage->dx, controlPointImage->dy, controlPointImage->dz);	
            printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
        }

		float maxStepSize = (targetImageCopy->dx > targetImageCopy->dy) ? targetImageCopy->dx : targetImageCopy->dy;
		maxStepSize = (targetImageCopy->dz > maxStepSize) ? targetImageCopy->dz : maxStepSize;

		float currentSize = maxStepSize;
		float smallestSize = maxStepSize / 100.0f;

        /* the gradient images are allocated */
        nifti_image *resultGradientImage = NULL;
        nifti_image *voxelNMIGradientImage = NULL;
        nifti_image *nodeNMIGradientImage = NULL;
        /* Conjugate gradient */
        PRECISION_TYPE *conjugateG = NULL;
        PRECISION_TYPE *conjugateH = NULL;
        /* joint histogram related variables */
        double *probaJointHistogram = (double *) malloc(nBins * (nBins+2) * sizeof(double));
        double *logJointHistogram = (double *) malloc(nBins * (nBins+2) * sizeof(double));
        double *entropies = (double *) malloc(4*sizeof(double));

        PRECISION_TYPE *bestControlPointPosition=NULL;

		resultGradientImage = nifti_copy_nim_info(positionFieldImage);
		if (sizeof(PRECISION_TYPE) == 4)
		    resultGradientImage->datatype = NIFTI_TYPE_FLOAT32;
		else
		    resultGradientImage->datatype = NIFTI_TYPE_FLOAT64;
		resultGradientImage->nbyper = sizeof(PRECISION_TYPE);
		resultGradientImage->data = (void *) calloc(resultGradientImage->nvox, resultGradientImage->nbyper);
		voxelNMIGradientImage = nifti_copy_nim_info(positionFieldImage);
		if (sizeof(PRECISION_TYPE) == 4)
		    voxelNMIGradientImage->datatype = NIFTI_TYPE_FLOAT32;
		else
		    voxelNMIGradientImage->datatype = NIFTI_TYPE_FLOAT64;
		voxelNMIGradientImage->nbyper = sizeof(PRECISION_TYPE);
		voxelNMIGradientImage->data = (void *) calloc(voxelNMIGradientImage->nvox, voxelNMIGradientImage->nbyper);
		nodeNMIGradientImage = nifti_copy_nim_info(controlPointImage);
		if (sizeof(PRECISION_TYPE) == 4)
		    nodeNMIGradientImage->datatype = NIFTI_TYPE_FLOAT32;
		else
		    nodeNMIGradientImage->datatype = NIFTI_TYPE_FLOAT64;
		nodeNMIGradientImage->nbyper = sizeof(PRECISION_TYPE);
		nodeNMIGradientImage->data = (void *) calloc(nodeNMIGradientImage->nvox, nodeNMIGradientImage->nbyper);
		
		conjugateG = (PRECISION_TYPE *)calloc(nodeNMIGradientImage->nvox, sizeof(PRECISION_TYPE));
		conjugateH = (PRECISION_TYPE *)calloc(nodeNMIGradientImage->nvox, sizeof(PRECISION_TYPE));
		
		bestControlPointPosition = (PRECISION_TYPE *) malloc(controlPointImage->nvox * sizeof(PRECISION_TYPE));

		memcpy(bestControlPointPosition, controlPointImage->data, controlPointImage->nvox*controlPointImage->nbyper);

		int smoothingRadius[3];
		smoothingRadius[0] = (int) floor(2.0*controlPointImage->dx/targetImageCopy->dx);
		smoothingRadius[1] = (int) floor(2.0*controlPointImage->dy/targetImageCopy->dy);
		smoothingRadius[2] = (int) floor(2.0*controlPointImage->dz/targetImageCopy->dz);

		int iteration = 0;

        while (iteration < maxIterations && currentSize > smallestSize)
        {
            double currentValue = 0.0;
            double currentWBE = 0.0f;
            double currentWJac = 0.0f;
            
            /* The Jacobian-based penalty term is first computed */
            if (jacobianWeight > 0)
                currentWJac = jacobianWeight * reg_bspline_jacobian<PRECISION_TYPE>(controlPointImage, targetImageCopy, true);

            /* The control point grid is corrected if necessary */
            int initialNegCorrection=0;
            while (currentWJac != currentWJac && initialNegCorrection < FOLDING_CORRECTION_STEP && iteration == 0)
            {
                currentWJac = jacobianWeight * reg_bspline_correctFolding<PRECISION_TYPE>(controlPointImage, targetImageCopy, true);
                initialNegCorrection++;
                if (currentWJac != currentWJac)
                {
                    if (verbose)
                        printf("*** Initial folding correction [%i/%i] ***\n", initialNegCorrection, FOLDING_CORRECTION_STEP);
                }
                else
                {
                    memcpy(bestControlPointPosition, controlPointImage->data, controlPointImage->nvox * controlPointImage->nbyper);
                    if (verbose)
                        printf(">>> Initial Jacobian based penalty term value = %g\n", currentWJac);
                }
            }
            if (currentWJac != currentWJac)
            {
                fprintf(stderr, "[WARNING] The initial folding correction failed.\n");
                fprintf(stderr, "[WARNING] You might want to increase the penalty term weights.\n");
                {
                    memcpy(controlPointImage->data, bestControlPointPosition, controlPointImage->nvox*controlPointImage->nbyper);
                }
            }

            /* The bending-energy penalty term is computed */
            if (bendingEnergyWeight > 0)
                currentWBE = bendingEnergyWeight * reg_bspline_bendingEnergy<PRECISION_TYPE>(controlPointImage, targetImageCopy, true);

            /* The source image is resampled and the metric assessed */

            /* generate the position field */
            reg_bspline<PRECISION_TYPE>(controlPointImage, targetImageCopy, positionFieldImage, targetMask, 0);
            
            /* Resample the source image */
            reg_resampleSourceImage<PRECISION_TYPE>(targetImageCopy, sourceImageCopy, resultImage, positionFieldImage, targetMask, 1, sourcePaddingValue);
            
            reg_getEntropies<double>(targetImageCopy, resultImage, JH_PW_APPROX, nBins, probaJointHistogram, logJointHistogram, entropies, targetMask);
			currentValue = (1.0 - bendingEnergyWeight - jacobianWeight) * (entropies[0] + entropies[1]) / entropies[2];

            iteration++;

            double bestWBE = currentWBE;
            double bestWJac = currentWJac;
            double bestValue = currentValue - bestWBE - currentWJac;
            if (iteration==1)
                printf("Initial objective function value = %g\n", bestValue);

            float maxLength;

			/* The NMI Gradient is calculated */
			reg_getSourceImageGradient<PRECISION_TYPE>(targetImageCopy, sourceImageCopy, resultGradientImage, positionFieldImage, targetMask, 1);
			reg_getVoxelBasedNMIGradientUsingPW<double>(targetImageCopy, resultImage, JH_PW_APPROX, resultGradientImage, nBins, logJointHistogram, entropies, voxelNMIGradientImage, targetMask);
            
            reg_smoothImageForCubicSpline<PRECISION_TYPE>(voxelNMIGradientImage, smoothingRadius);
            reg_voxelCentric2NodeCentric(nodeNMIGradientImage, voxelNMIGradientImage, 1.0f - bendingEnergyWeight - jacobianWeight);

            /* The NMI gradient is converted from voxel space to real space */
            if (twoDimRegistration)
            {
                PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
                PRECISION_TYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx * controlPointImage->ny * controlPointImage->nz];
                PRECISION_TYPE newGradientValueX, newGradientValueY;
                for (int i=0; i < controlPointImage->nx * controlPointImage->ny * controlPointImage->nz; i++)
                {
	                newGradientValueX = *gradientValuesX * sourceMatrix_xyz->m[0][0] + *gradientValuesY * sourceMatrix_xyz->m[0][1];
	                newGradientValueY = *gradientValuesX * sourceMatrix_xyz->m[1][0] + *gradientValuesY * sourceMatrix_xyz->m[1][1];
	                *gradientValuesX++ = newGradientValueX;
	                *gradientValuesY++ = newGradientValueY;
                }
            }
            else
            {
                PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
                PRECISION_TYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx * controlPointImage->ny * controlPointImage->nz];
                PRECISION_TYPE *gradientValuesZ = &gradientValuesY[controlPointImage->nx * controlPointImage->ny * controlPointImage->nz];
                PRECISION_TYPE newGradientValueX, newGradientValueY, newGradientValueZ;
                for (int i=0; i < controlPointImage->nx * controlPointImage->ny * controlPointImage->nz; i++)
                {
                    newGradientValueX = *gradientValuesX * reorient_NMI_gradient.m[0][0] + *gradientValuesY * reorient_NMI_gradient.m[0][1] + *gradientValuesZ * reorient_NMI_gradient.m[0][2];
                    newGradientValueY = *gradientValuesX * reorient_NMI_gradient.m[1][0] + *gradientValuesY * reorient_NMI_gradient.m[1][1] + *gradientValuesZ * reorient_NMI_gradient.m[1][2];
                    newGradientValueZ = *gradientValuesX * reorient_NMI_gradient.m[2][0] + *gradientValuesY * reorient_NMI_gradient.m[2][1] + *gradientValuesZ * reorient_NMI_gradient.m[2][2];
                    *gradientValuesX++ = newGradientValueX;
                    *gradientValuesY++ = newGradientValueY;
                    *gradientValuesZ++ = newGradientValueZ;
                }
            }
            
            /* The other gradients are calculated */
            if (bendingEnergyWeight > 0)
				reg_bspline_bendingEnergyGradient<PRECISION_TYPE>(controlPointImage, targetImageCopy, nodeNMIGradientImage, bendingEnergyWeight);
				
			if(jacobianWeight > 0)
                reg_bspline_jacobianDeterminantGradient<PRECISION_TYPE>(controlPointImage, targetImageCopy, nodeNMIGradientImage, jacobianWeight, true);

			/* The conjugate gradient is computed */
			if (iteration==1)
			{
				// first conjugate gradient iteration
				if (twoDimRegistration)
				{
					PRECISION_TYPE *conjGPtrX = &conjugateG[0];
					PRECISION_TYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
					PRECISION_TYPE *conjHPtrX = &conjugateH[0];
					PRECISION_TYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
					PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
					PRECISION_TYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny];
					for (int i=0; i < nodeNMIGradientImage->nx * nodeNMIGradientImage->ny; i++)
					{
						*conjHPtrX++ = *conjGPtrX++ = - *gradientValuesX++;
						*conjHPtrY++ = *conjGPtrY++ = - *gradientValuesY++;
					}
				}
				else
				{
					PRECISION_TYPE *conjGPtrX = &conjugateG[0];
					PRECISION_TYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *conjGPtrZ = &conjGPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *conjHPtrX = &conjugateH[0];
					PRECISION_TYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *conjHPtrZ = &conjHPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
					PRECISION_TYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
					PRECISION_TYPE *gradientValuesZ = &gradientValuesY[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
					for(int i=0; i < nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz; i++)
					{
						*conjHPtrX++ = *conjGPtrX++ = - *gradientValuesX++;
						*conjHPtrY++ = *conjGPtrY++ = - *gradientValuesY++;
						*conjHPtrZ++ = *conjGPtrZ++ = - *gradientValuesZ++;
					}
				}
			}
			else
			{
				double dgg = 0.0, gg = 0.0;
				if (twoDimRegistration)
				{
					PRECISION_TYPE *conjGPtrX = &conjugateG[0];
					PRECISION_TYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
					PRECISION_TYPE *conjHPtrX = &conjugateH[0];
					PRECISION_TYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
					PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
					PRECISION_TYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny];
					for(int i=0; i < nodeNMIGradientImage->nx * nodeNMIGradientImage->ny; i++)
					{
						gg += conjHPtrX[i] * conjGPtrX[i];
						gg += conjHPtrY[i] * conjGPtrY[i];
						dgg += (gradientValuesX[i] + conjGPtrX[i]) * gradientValuesX[i];
						dgg += (gradientValuesY[i] + conjGPtrY[i]) * gradientValuesY[i];
					}
					double gam = dgg/gg;
					for (int i = 0; i < nodeNMIGradientImage->nx * nodeNMIGradientImage->ny; i++)
					{
						conjGPtrX[i] = -gradientValuesX[i];
						conjGPtrY[i] = -gradientValuesY[i];
						conjHPtrX[i] = (float)(conjGPtrX[i] + gam * conjHPtrX[i]);
						conjHPtrY[i] = (float)(conjGPtrY[i] + gam * conjHPtrY[i]);
						gradientValuesX[i] = -conjHPtrX[i];
						gradientValuesY[i] = -conjHPtrY[i];
					}
				}
				else
				{
					PRECISION_TYPE *conjGPtrX = &conjugateG[0];
					PRECISION_TYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *conjGPtrZ = &conjGPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *conjHPtrX = &conjugateH[0];
					PRECISION_TYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *conjHPtrZ = &conjHPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
					PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
					PRECISION_TYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
					PRECISION_TYPE *gradientValuesZ = &gradientValuesY[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
					for (int i = 0; i < nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz; i++)
					{
						gg += conjHPtrX[i] * conjGPtrX[i];
						gg += conjHPtrY[i] * conjGPtrY[i];
						gg += conjHPtrZ[i] * conjGPtrZ[i];
						dgg += (gradientValuesX[i] + conjGPtrX[i]) * gradientValuesX[i];
						dgg += (gradientValuesY[i] + conjGPtrY[i]) * gradientValuesY[i];
						dgg += (gradientValuesZ[i] + conjGPtrZ[i]) * gradientValuesZ[i];
					}
					double gam = dgg/gg;
					for (int i = 0; i < nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz; i++)
					{
						conjGPtrX[i] = -gradientValuesX[i];
						conjGPtrY[i] = -gradientValuesY[i];
						conjGPtrZ[i] = -gradientValuesZ[i];
						conjHPtrX[i] = (float)(conjGPtrX[i] + gam * conjHPtrX[i]);
						conjHPtrY[i] = (float)(conjGPtrY[i] + gam * conjHPtrY[i]);
						conjHPtrZ[i] = (float)(conjGPtrZ[i] + gam * conjHPtrZ[i]);
						gradientValuesX[i] = -conjHPtrX[i];
						gradientValuesY[i] = -conjHPtrY[i];
						gradientValuesZ[i] = -conjHPtrZ[i];
					}
				}
			}
			maxLength = reg_getMaximalLength<PRECISION_TYPE>(nodeNMIGradientImage);

			/* The gradient is applied to the control point positions */
			if (maxLength==0)
			{
				printf("No Gradient ... exit\n");
				break;	
			}

			/* ** LINE ASCENT ** */
			int lineIteration = 0;
    		currentSize = maxStepSize;
			float addedStep = 0.0f;

            while (currentSize > smallestSize && lineIteration < 12)
            {

				float currentLength = -currentSize/maxLength;

                // the control point positions are updated using addition
                if (twoDimRegistration)
                {
                    PRECISION_TYPE *controlPointValuesX = NULL;
                    PRECISION_TYPE *controlPointValuesY = NULL;
                    controlPointValuesX = static_cast<PRECISION_TYPE *>(controlPointImage->data);
                    controlPointValuesY = &controlPointValuesX[controlPointImage->nx*controlPointImage->ny];
                    PRECISION_TYPE *bestControlPointValuesX = &bestControlPointPosition[0];
                    PRECISION_TYPE *bestControlPointValuesY = &bestControlPointValuesX[controlPointImage->nx*controlPointImage->ny];
                    PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
                    PRECISION_TYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx*controlPointImage->ny];
                    for (int i = 0; i < controlPointImage->nx * controlPointImage->ny; i++)
                    {
                        *controlPointValuesX++ = *bestControlPointValuesX++ + currentLength * *gradientValuesX++;
                        *controlPointValuesY++ = *bestControlPointValuesY++ + currentLength * *gradientValuesY++;
                    }
                }
                else
                {
                    PRECISION_TYPE *controlPointValuesX = NULL;
                    PRECISION_TYPE *controlPointValuesY = NULL;
                    PRECISION_TYPE *controlPointValuesZ = NULL;
                    controlPointValuesX = static_cast<PRECISION_TYPE *>(controlPointImage->data);
                    controlPointValuesY = &controlPointValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                    controlPointValuesZ = &controlPointValuesY[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                    PRECISION_TYPE *bestControlPointValuesX = &bestControlPointPosition[0];
                    PRECISION_TYPE *bestControlPointValuesY = &bestControlPointValuesX[controlPointImage->nx * controlPointImage->ny * controlPointImage->nz];
                    PRECISION_TYPE *bestControlPointValuesZ = &bestControlPointValuesY[controlPointImage->nx * controlPointImage->ny * controlPointImage->nz];
                    PRECISION_TYPE *gradientValuesX = static_cast<PRECISION_TYPE *>(nodeNMIGradientImage->data);
                    PRECISION_TYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                    PRECISION_TYPE *gradientValuesZ = &gradientValuesY[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                    for(int i=0; i<controlPointImage->nx*controlPointImage->ny*controlPointImage->nz;i++)
                    {
                        *controlPointValuesX++ = *bestControlPointValuesX++ + currentLength * *gradientValuesX++;
                        *controlPointValuesY++ = *bestControlPointValuesY++ + currentLength * *gradientValuesY++;
                        *controlPointValuesZ++ = *bestControlPointValuesZ++ + currentLength * *gradientValuesZ++;
                    }
                }

                /* The Jacobian-based penalty term is computed */
                if (jacobianWeight > 0)
                    currentWJac = jacobianWeight * reg_bspline_jacobian<PRECISION_TYPE>(controlPointImage, targetImageCopy, true);

                int negCorrection = 0;
                while (currentWJac != currentWJac && negCorrection < 5)
                {
                    printf("*");
                    currentWJac = jacobianWeight * reg_bspline_correctFolding<PRECISION_TYPE>(controlPointImage, targetImageCopy, true);
                    negCorrection++;
                }

                /* The bending energy is computed */
                if (bendingEnergyWeight > 0)
                    currentWBE = bendingEnergyWeight * reg_bspline_bendingEnergy<PRECISION_TYPE>(controlPointImage, targetImageCopy, true);

                /* The source image is resampled and the metric evaluated */
				reg_bspline<PRECISION_TYPE>(controlPointImage, targetImageCopy, positionFieldImage, targetMask, 0);

				/* Resample the source image */
				reg_resampleSourceImage<PRECISION_TYPE>(targetImageCopy, sourceImageCopy, resultImage, positionFieldImage, NULL, 1, sourcePaddingValue);
				
				/* Computation of the Metric Value */
				reg_getEntropies<double>(targetImageCopy, resultImage, JH_PW_APPROX, nBins, probaJointHistogram, logJointHistogram, entropies, targetMask);
                currentValue = (1.0 - bendingEnergyWeight - jacobianWeight) * (entropies[0] + entropies[1]) / entropies[2];

                currentValue -= currentWBE + currentWJac;

                iteration++;
                lineIteration++;

				/* The current deformation field is kept if it was the best so far */
				if (currentValue > bestValue)
				{
					memcpy(bestControlPointPosition, controlPointImage->data, controlPointImage->nvox * controlPointImage->nbyper);
					bestValue = currentValue;
					bestWBE = currentWBE;
					bestWJac = currentWJac;
					addedStep += currentSize;
					currentSize *= 1.1f;
					currentSize = (currentSize < maxStepSize) ? currentSize : maxStepSize;
				}
				else
					currentSize *= 0.5;
			}
			
			memcpy(controlPointImage->data, bestControlPointPosition, controlPointImage->nvox * controlPointImage->nbyper);
			currentSize = addedStep;
			printf("[%i] Objective function value=%g | max. added disp. = %g mm", iteration, bestValue, addedStep);
			if (bendingEnergyWeight > 0)
			    printf(" | wBE=%g", bestWBE);
			if (jacobianWeight > 0)
			    printf(" | wJacLog=%g", bestWJac);
			printf("\n");
		}  // end of interation loop
		
		completedIterations[level] = iteration;

        /* The deformation model is unfolded if the Jacobian determinant-based penalty term is used */
        if (jacobianWeight > 0)
        {
            int finalNegCorrection = 0;
            PRECISION_TYPE finalWJac = 0;
            do
            {
                finalWJac = jacobianWeight * reg_bspline_correctFolding<PRECISION_TYPE>(controlPointImage, targetImageCopy, false); 
                finalNegCorrection++;
                if (finalWJac != finalWJac)
                    printf( "*** Final folding correction [%i/%i] ***\n", finalNegCorrection, FOLDING_CORRECTION_STEP);
                else
                    printf(">>> Final Jacobian based penalty term value = %g\n", finalWJac);
            }
            while (finalWJac != finalWJac && finalNegCorrection < FOLDING_CORRECTION_STEP);

            if(finalWJac != finalWJac)
            {
                fprintf(stderr, "[WARNING] The final folding correction failed.\n");
                fprintf(stderr, "[WARNING] You might want to increase the penalty term weights.\n");
                {
                    memcpy(controlPointImage->data, bestControlPointPosition, controlPointImage->nvox * controlPointImage->nbyper);
                }
            }
        }

        free(targetMask);
		free(entropies);
		free(probaJointHistogram);
		free(logJointHistogram);
		
		nifti_image_free(resultGradientImage);
		nifti_image_free(voxelNMIGradientImage);
		nifti_image_free(nodeNMIGradientImage);
		
		free(conjugateG);
		free(conjugateH);
        free(bestControlPointPosition);
        
        nifti_image_free(resultImage);
        nifti_image_free(targetImageCopy);
        nifti_image_free(sourceImageCopy);
        
        if (level < (nLevels - 1))
            nifti_image_free(positionFieldImage);
        
        if (verbose)
    		printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
	} // end of level loop
	
	if (nLevels == 0)
        positionFieldImage = create_position_field(targetImage, twoDimRegistration);
	
	reg_bspline<PRECISION_TYPE>(controlPointImage, targetImage, positionFieldImage, NULL, 0);
	
	// The source data type is changed for precision if requested
    if (finalPrecision == INTERP_PREC_DOUBLE)
        reg_changeDatatype<double>(sourceImage);
    
    // The result image is resampled using a cubic spline interpolation
    resultImage = nifti_copy_nim_info(targetImage);
    resultImage->cal_min = sourceImage->cal_min;
    resultImage->cal_max = sourceImage->cal_max;
    resultImage->scl_slope = sourceImage->scl_slope;
    resultImage->scl_inter = sourceImage->scl_inter;
    resultImage->datatype = sourceImage->datatype;
    resultImage->nbyper = sourceImage->nbyper;
    resultImage->data = calloc(resultImage->nvox, resultImage->nbyper);
    reg_resampleSourceImage<PRECISION_TYPE>(targetImage, sourceImage, resultImage, positionFieldImage, NULL, finalInterpolation, sourcePaddingValue);
    
    nifti_image_free(positionFieldImage);
    
    f3d_result result;
    result.image = resultImage;
    result.controlPoints = controlPointImage;
    result.completedIterations = completedIterations;
    
    return result;
}
