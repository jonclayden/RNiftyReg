// Registration types (degrees of freedom)
#define RIGID 0
#define AFFINE 1

// Preceision levels for final image
#define SOURCE 0
#define SINGLE 1
#define DOUBLE 2

// Working precision for the source and target images: float and double are valid for the NiftyReg code
#define PRECISION_TYPE double

// Convergence criterion
#define CONVERGENCE_EPS 0.00001

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_blockMatching.h"
#include "_reg_tools.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "niftyreg.h"

extern "C"
SEXP reg_aladin (SEXP source, SEXP target, SEXP type, SEXP finalPrecision, SEXP nLevels, SEXP maxIterations, SEXP useBlockPercentage, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP verbose)
{
    int i, j;
    SEXP returnValue, data;
    
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
    const char *precisionString = CHAR(STRING_ELT(finalPrecision,0));
    int precisionType = SOURCE;
    
    if (strcmp(precisionString,"single") == 0)
        precisionType = SINGLE;
    else if (strcmp(precisionString,"double") == 0)
        precisionType = DOUBLE;
    
    aladin_result result = do_reg_aladin(sourceImage, targetImage, regType, precisionType, *(INTEGER(nLevels)), *(INTEGER(maxIterations)), *(INTEGER(useBlockPercentage)), *(INTEGER(finalInterpolation)), targetMaskImage, affineTransformation, (*(INTEGER(verbose)) == 1));
    
    PROTECT(returnValue = NEW_LIST(2));
    
    // Integer-valued data went in, and precision must be "source"
    if (result.image->datatype == DT_INT32)
    {
        PROTECT(data = NEW_INTEGER((R_len_t) result.image->nvox));
        for (size_t i = 0; i < result.image->nvox; i++)
            INTEGER(data)[i] = ((int *) result.image->data)[i];
    }
    else
    {
        PROTECT(data = NEW_NUMERIC((R_len_t) result.image->nvox));
        if (precisionType == SINGLE && *(INTEGER(finalInterpolation)) != 0)
        {
            for (size_t i = 0; i < result.image->nvox; i++)
                REAL(data)[i] = (double) ((float *) result.image->data)[i];
        }
        else
        {
            for (size_t i = 0; i < result.image->nvox; i++)
                REAL(data)[i] = ((double *) result.image->data)[i];
        }
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
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    if (targetMaskImage != NULL)
        nifti_image_free(targetMaskImage);
    nifti_image_free(result.image);
    free(result.affine);
    
    UNPROTECT(affineProvided ? 2 : 3);
    
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
        fprintf(stderr, "** ERROR: Data type %d is not supported\n", header.datatype);
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
            printf("Current level %i / %i\n", level+1, nLevels);
            printf("Target image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n", targetImageCopy->nx, targetImageCopy->ny, targetImageCopy->nz, targetImageCopy->dx, targetImageCopy->dy, targetImageCopy->dz);
            printf("Source image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n", sourceImageCopy->nx, sourceImageCopy->ny, sourceImageCopy->nz, sourceImageCopy->dx, sourceImageCopy->dy, sourceImageCopy->dz);
            if (twoDimRegistration)
                printf("Block size = [4 4 1]\n");
            else
                printf("Block size = [4 4 4]\n");
            printf("Block number = [%i %i %i]\n", blockMatchingParams.blockNumber[0], blockMatchingParams.blockNumber[1], blockMatchingParams.blockNumber[2]);
            printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
            reg_mat44_disp(affineTransformation, (char *) "Initial affine transformation");
        }

        int nLoops = ((type==AFFINE && level==0) ? 2 : 1);
        int currentType, iteration = 0;
        mat44 updateAffineMatrix;
        
        for (i = 0; i < nLoops; i++)
        {
            currentType = ((i==0 && nLoops==2) ? RIGID : type);
            
            // Twice as many iterations are performed during the first level
            while (iteration < (level==0 ? maxIterations : 2*maxIterations))
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

        free(targetMask);
        nifti_image_free(resultImage);
        nifti_image_free(targetImageCopy);
        nifti_image_free(sourceImageCopy);
        
        if (level < (nLevels - 1))
            nifti_image_free(positionFieldImage);
        
        if (verbose)
        {
            reg_mat44_disp(affineTransformation, (char *)"Final affine transformation");
            printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
        }
    }
    
    if (nLevels == 0)
        positionFieldImage = create_position_field(targetImage, twoDimRegistration);
    
    // The corresponding deformation field is evaluated and saved
    reg_affine_positionField(affineTransformation, targetImage, positionFieldImage);
    
    // The source data type is changed to ensure precision unless using nearest neighbour interpolation is to be used
    if (finalPrecision == SINGLE && finalInterpolation != 0)
        reg_changeDatatype<float>(sourceImage);
    else if (finalPrecision == DOUBLE && finalInterpolation != 0)
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
    
    return result;
}
