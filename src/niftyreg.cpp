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

f3d_result do_reg_f3d (nifti_image *sourceImage, nifti_image *targetImage, int type, int finalPrecision, int nLevels, int maxIterations, int useBlockPercentage, int finalInterpolation, nifti_image *targetMaskImage, nifti_image *controlPointImage, mat44 *affineTransformation, int nBins, bool verbose)
{
    bool twoDimRegistration = (sourceImage->nz == 1 || targetImage->nz == 1);
    
    // This is due to the extrapolation of the joint histogram using the Parzen window
    nBins += 4;
    
    // param->binning += 4; //This is due to the extrapolation of the joint histogram using the Parzen window
    // if(param->spacing[0]<0) param->spacing[0] *=
    //     -1.0f * targetHeader->dx * powf(2.0f, (float)(param->levelNumber-param->level2Perform));
    // if(param->spacing[1]<1) param->spacing[1] *=
    //     -1.0f * targetHeader->dy * powf(2.0f, (float)(param->levelNumber-param->level2Perform));
    // if(param->spacing[2]<2) param->spacing[2] *=
    //     -1.0f * targetHeader->dz * powf(2.0f, (float)(param->levelNumber-param->level2Perform));
    // param->sourcePaddingValue = std::numeric_limits<float>::quiet_NaN();
    // nifti_image *controlPointImage=NULL;
    
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
            if (controlPointImage != NULL)
            {
                // Allocate the control point image
                int dim_cpp[8];
                float gridSpacing[3];
                dim_cpp[0] = 5;
                gridSpacing[0] = param->spacing[0] * powf(2.0f, (float)(param->level2Perform-1));
                dim_cpp[1] = (int) floor(targetImageCopy->nx*targetImageCopy->dx/gridSpacing[0]) + 5;
                gridSpacing[1] = param->spacing[1] * powf(2.0f, (float)(param->level2Perform-1));
                dim_cpp[2] = (int) floor(targetImageCopy->ny*targetImageCopy->dy/gridSpacing[1]) + 5;
                if (twoDimRegistration)
                {
                    gridSpacing[2] = 1.0f;
                	dim_cpp[3]=1;
                	dim_cpp[5]=2;
                }
                else
                {
                    gridSpacing[2] = param->spacing[2] * powf(2.0f, (float)(param->level2Perform-1));
                    dim_cpp[3]=(int)floor(targetImageCopy->nz*targetImageCopy->dz/gridSpacing[2])+5;
                	dim_cpp[5]=3;
                }
                dim_cpp[4] = dim_cpp[6] = dim_cpp[7] = 1;
                
                if (sizeof(PrecisionTYPE) == 4)
                    controlPointImage = nifti_make_new_nim(dim_cpp, NIFTI_TYPE_FLOAT32, true);
                else
                    controlPointImage = nifti_make_new_nim(dim_cpp, NIFTI_TYPE_FLOAT64, true);
                
                controlPointImage->cal_min = 0;
                controlPointImage->cal_max = 0;
                controlPointImage->pixdim[0] = 1.0f;
                controlPointImage->pixdim[1] = controlPointImage->dx=gridSpacing[0];
                controlPointImage->pixdim[2] = controlPointImage->dy=gridSpacing[1];
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
            reg_bspline_refineControlPointGrid(targetImage, controlPointImage);

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
            controlPointImage->qform_code=1;
        controlPointImage->qto_xyz.m[0][3] = controlPointImage->qoffset_x = originReal[0];
        controlPointImage->qto_xyz.m[1][3] = controlPointImage->qoffset_y = originReal[1];
        controlPointImage->qto_xyz.m[2][3] = controlPointImage->qoffset_z = originReal[2];

        controlPointImage->qto_ijk = nifti_mat44_inverse(controlPointImage->qto_xyz);

        if (controlPointImage->sform_code > 0)
        {
			nifti_mat44_to_quatern(targetImage->sto_xyz, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac);

			controlPointImage->sto_xyz = nifti_quatern_to_mat44(qb, qc, qd, qx, qy, qz, controlPointImage->dx, controlPointImage->dy, controlPointImage->dz, qfac);

			// Origin is shifted from 1 control point in the sform
			originIndex[0] = -1.0f;
			originIndex[1] = -1.0f;
			originIndex[2] = 0.0f;
			if(targetImage->nz>1) originIndex[2] = -1.0f;
			reg_mat44_mul(&(controlPointImage->sto_xyz), originIndex, originReal);
			controlPointImage->sto_xyz.m[0][3] = originReal[0];
			controlPointImage->sto_xyz.m[1][3] = originReal[1];
			controlPointImage->sto_xyz.m[2][3] = originReal[2];

            controlPointImage->sto_ijk = nifti_mat44_inverse(controlPointImage->sto_xyz);
        }

        if (level==0 && !flag->inputCPPFlag)
        {
            // The control point position image is initialised with the affine transformation
            if (reg_bspline_initialiseControlPointGridWithAffine(affineTransformation, controlPointImage))
                return 1;
            free(affineTransformation);
        }

        mat44 *cppMatrix_xyz;
        mat44 *targetMatrix_ijk;
        mat44 *sourceMatrix_xyz;
        if (controlPointImage->sform_code)
            cppMatrix_xyz = &(controlPointImage->sto_xyz);
        else
            cppMatrix_xyz = &(controlPointImage->qto_xyz);
        if (targetImage->sform_code)
            targetMatrix_ijk = &(targetImage->sto_ijk);
        else
            targetMatrix_ijk = &(targetImage->qto_ijk);
        if (sourceImage->sform_code)
            sourceMatrix_xyz = &(sourceImage->sto_xyz);
        else
            sourceMatrix_xyz = &(sourceImage->qto_xyz);

        mat33 reorient_NMI_gradient;
        for(unsigned int i=0; i<3; i++)
        {
            for(unsigned int  j=0; j<3; j++)
                reorient_NMI_gradient.m[i][j] = sourceMatrix_xyz->m[i][j];
        }

        /* allocate the deformation Field image */
        nifti_image *positionFieldImage = nifti_copy_nim_info(targetImage);
        positionFieldImage->dim[0] = positionFieldImage->ndim = 5;
        positionFieldImage->dim[1] = positionFieldImage->nx = targetImage->nx;
        positionFieldImage->dim[2] = positionFieldImage->ny = targetImage->ny;
        positionFieldImage->dim[3] = positionFieldImage->nz = targetImage->nz;
        positionFieldImage->dim[4] = positionFieldImage->nt = 1;
        positionFieldImage->pixdim[4] = positionFieldImage->dt = 1.0;
        if (twoDimRegistration)
            positionFieldImage->dim[5] = positionFieldImage->nu = 2;
        else
            positionFieldImage->dim[5] = positionFieldImage->nu = 3;
        positionFieldImage->pixdim[5] = positionFieldImage->du = 1.0;
        positionFieldImage->dim[6] = positionFieldImage->nv = 1;
        positionFieldImage->pixdim[6] = positionFieldImage->dv = 1.0;
        positionFieldImage->dim[7] = positionFieldImage->nw = 1;
        positionFieldImage->pixdim[7] = positionFieldImage->dw = 1.0;
        positionFieldImage->nvox = positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
        if (sizeof(PrecisionTYPE)==4)
            positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
        else
            positionFieldImage->datatype = NIFTI_TYPE_FLOAT64;
        positionFieldImage->nbyper = sizeof(PrecisionTYPE);
        positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);

        /* allocate the result image */
        nifti_image *resultImage = nifti_copy_nim_info(targetImage);
        resultImage->datatype = sourceImage->datatype;
        resultImage->nbyper = sourceImage->nbyper;
        resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);

        if (verbose)
        {
    		printf("Current level %i / %i\n", level+1, param->levelNumber);
    		printf("Target image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n",
    		       targetImage->nx, targetImage->ny, targetImage->nz, targetImage->dx, targetImage->dy, targetImage->dz);
    		printf("Source image size: \t%ix%ix%i voxels\t%gx%gx%g mm\n",
    		       sourceImage->nx, sourceImage->ny, sourceImage->nz, sourceImage->dx, sourceImage->dy, sourceImage->dz);
    		printf("Control point position image name: %s\n",param->outputCPPName);
    		printf("\t%ix%ix%i control points (%i DoF)\n",controlPointImage->nx,controlPointImage->ny,controlPointImage->nz,(int)controlPointImage->nvox);
    		printf("\t%gx%gx%g mm\n",controlPointImage->dx,controlPointImage->dy,controlPointImage->dz);	
            printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
        }

		float maxStepSize = (targetImage->dx>targetImage->dy)?targetImage->dx:targetImage->dy;
		maxStepSize = (targetImage->dz>maxStepSize)?targetImage->dz:maxStepSize;

		float currentSize = maxStepSize;
		float smallestSize = maxStepSize / 100.0f;

        if(flag->backgroundIndexFlag){
            int index[3];
            index[0]=param->backgroundIndex[0];
            index[1]=param->backgroundIndex[1];
            index[2]=param->backgroundIndex[2];
            if(flag->pyramidFlag){
                for(int l=level; l<param->levelNumber-1; l++){
                    index[0] /= 2;
                    index[1] /= 2;
                    index[2] /= 2;
                }
            }
        }

        /* the gradient images are allocated */
        nifti_image *resultGradientImage=NULL;
        nifti_image *voxelNMIGradientImage=NULL;
        nifti_image *nodeNMIGradientImage=NULL;
        /* Conjugate gradient */
        PrecisionTYPE *conjugateG=NULL;
        PrecisionTYPE *conjugateH=NULL;
        /* joint histogram related variables */
        double *probaJointHistogram = (double *)malloc(param->binning*(param->binning+2)*sizeof(double));
        double *logJointHistogram = (double *)malloc(param->binning*(param->binning+2)*sizeof(double));
        double *entropies = (double *)malloc(4*sizeof(double));

        PrecisionTYPE *bestControlPointPosition=NULL;

			resultGradientImage = nifti_copy_nim_info(positionFieldImage);
			if(sizeof(PrecisionTYPE)==4) resultGradientImage->datatype = NIFTI_TYPE_FLOAT32;
			else resultGradientImage->datatype = NIFTI_TYPE_FLOAT64;
			resultGradientImage->nbyper = sizeof(PrecisionTYPE);
			resultGradientImage->data = (void *)calloc(resultGradientImage->nvox, resultGradientImage->nbyper);
			voxelNMIGradientImage = nifti_copy_nim_info(positionFieldImage);
			if(sizeof(PrecisionTYPE)==4) voxelNMIGradientImage->datatype = NIFTI_TYPE_FLOAT32;
			else voxelNMIGradientImage->datatype = NIFTI_TYPE_FLOAT64;
			voxelNMIGradientImage->nbyper = sizeof(PrecisionTYPE);
			voxelNMIGradientImage->data = (void *)calloc(voxelNMIGradientImage->nvox, voxelNMIGradientImage->nbyper);
			nodeNMIGradientImage = nifti_copy_nim_info(controlPointImage);
			if(sizeof(PrecisionTYPE)==4) nodeNMIGradientImage->datatype = NIFTI_TYPE_FLOAT32;
			else nodeNMIGradientImage->datatype = NIFTI_TYPE_FLOAT64;
			nodeNMIGradientImage->nbyper = sizeof(PrecisionTYPE);
			nodeNMIGradientImage->data = (void *)calloc(nodeNMIGradientImage->nvox, nodeNMIGradientImage->nbyper);
			if(!flag->noConjugateGradient){
				conjugateG = (PrecisionTYPE *)calloc(nodeNMIGradientImage->nvox, sizeof(PrecisionTYPE));
				conjugateH = (PrecisionTYPE *)calloc(nodeNMIGradientImage->nvox, sizeof(PrecisionTYPE));
			}
			bestControlPointPosition = (PrecisionTYPE *)malloc(controlPointImage->nvox * sizeof(PrecisionTYPE));

			memcpy(bestControlPointPosition, controlPointImage->data, controlPointImage->nvox*controlPointImage->nbyper);

		int smoothingRadius[3];
		smoothingRadius[0] = (int)floor( 2.0*controlPointImage->dx/targetImage->dx );
		smoothingRadius[1] = (int)floor( 2.0*controlPointImage->dy/targetImage->dy );
		smoothingRadius[2] = (int)floor( 2.0*controlPointImage->dz/targetImage->dz );


		int iteration=0;

        while(iteration<param->maxIteration && currentSize>smallestSize){

            double currentValue=0.0;
            double currentWBE=0.0f;
            double currentWJac=0.0f;
            PrecisionTYPE SSDValue=0.0;

            /* The Jacobian-based penalty term is first computed */
                if(flag->jacobianWeightFlag && param->jacobianWeight>0){
                    currentWJac = param->jacobianWeight
                        * reg_bspline_jacobian<PrecisionTYPE>(controlPointImage, targetImage, flag->appJacobianFlag);
                }

            /* The control point grid is corrected if necessary */
            int initialNegCorrection=0;
            while(currentWJac!=currentWJac && initialNegCorrection<FOLDING_CORRECTION_STEP && iteration==0){
                    currentWJac = param->jacobianWeight*
                          reg_bspline_correctFolding<PrecisionTYPE>(controlPointImage,
                                                                    targetImage,
                                                                    flag->appJacobianFlag);
                initialNegCorrection++;
                if(currentWJac!=currentWJac){
                    printf("*** Initial folding correction [%i/%i] ***\n",
                       initialNegCorrection, FOLDING_CORRECTION_STEP);
                }
                else{
                        memcpy(bestControlPointPosition,controlPointImage->data,controlPointImage->nvox*controlPointImage->nbyper);
                    printf(">>> Initial Jacobian based penalty term value = %g\n", currentWJac);
                }
            }
            if(currentWJac!=currentWJac){
                fprintf(stderr, "[WARNING] The initial folding correction failed.\n");
                fprintf(stderr, "[WARNING] You might want to increase the penalty term weights.\n");
                {
                    memcpy(controlPointImage->data,bestControlPointPosition,controlPointImage->nvox*controlPointImage->nbyper);
                }
            }

            /* The bending-energy penalty term is computed */
                if(flag->bendingEnergyFlag && param->bendingEnergyWeight>0){
                    currentWBE = param->bendingEnergyWeight
                        * reg_bspline_bendingEnergy<PrecisionTYPE>(controlPointImage, targetImage, flag->appBendingEnergyFlag);
                }

            /* The source image is resampled and the metric assessed */

                /* generate the position field */
                reg_bspline<PrecisionTYPE>( controlPointImage,
                                            targetImage,
                                            positionFieldImage,
                                            targetMask,
                                            0);
                /* Resample the source image */
                reg_resampleSourceImage<PrecisionTYPE>(	targetImage,
                                                        sourceImage,
                                                        resultImage,
                                                        positionFieldImage,
                                                        targetMask,
                                                        1,
                                                        param->sourcePaddingValue);
            if(flag->useSSDFlag){
                SSDValue = reg_getSSD<PrecisionTYPE>(	targetImage,
                                                        resultImage);
                currentValue = -(1.0-param->bendingEnergyWeight-param->jacobianWeight) * log(SSDValue+1.0);
			}
			else{
                reg_getEntropies<double>(	targetImage,
								resultImage,
								JH_PW_APPROX,
								param->binning,
								probaJointHistogram,
								logJointHistogram,
								entropies,
								targetMask);
				currentValue = (1.0-param->bendingEnergyWeight-param->jacobianWeight)*(entropies[0]+entropies[1])/entropies[2];
			}



            iteration++;


            double bestWBE = currentWBE;
            double bestWJac = currentWJac;
            double bestValue = currentValue - bestWBE - currentWJac;
            if(iteration==1) printf("Initial objective function value = %g\n", bestValue);

            float maxLength;

				/* The NMI Gradient is calculated */
				reg_getSourceImageGradient<PrecisionTYPE>(	targetImage,
															sourceImage,
															resultGradientImage,
															positionFieldImage,
															targetMask,
															1);
				if(flag->useSSDFlag){
					reg_getVoxelBasedSSDGradient<PrecisionTYPE>(SSDValue,
																targetImage,
																resultImage,
																resultGradientImage,
																voxelNMIGradientImage);
				}
				else{
                    reg_getVoxelBasedNMIGradientUsingPW<double>(targetImage,
																resultImage,
																JH_PW_APPROX,
																resultGradientImage,
																param->binning,
																logJointHistogram,
																entropies,
																voxelNMIGradientImage,
																targetMask);
				}
                reg_smoothImageForCubicSpline<PrecisionTYPE>(voxelNMIGradientImage,smoothingRadius);
                reg_voxelCentric2NodeCentric(nodeNMIGradientImage,
											 voxelNMIGradientImage,
											 1.0f-param->bendingEnergyWeight-param->jacobianWeight);

                /* The NMI gradient is converted from voxel space to real space */
                if(flag->twoDimRegistration){
	                PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
	                PrecisionTYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
	                PrecisionTYPE newGradientValueX, newGradientValueY;
	                for(int i=0; i<controlPointImage->nx*controlPointImage->ny*controlPointImage->nz; i++){

		                newGradientValueX = 	*gradientValuesX * sourceMatrix_xyz->m[0][0] +
					                *gradientValuesY * sourceMatrix_xyz->m[0][1];
		                newGradientValueY = 	*gradientValuesX * sourceMatrix_xyz->m[1][0] +
					                *gradientValuesY * sourceMatrix_xyz->m[1][1];

		                *gradientValuesX++ = newGradientValueX;
		                *gradientValuesY++ = newGradientValueY;
	                }
                }
                else{
	                PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
	                PrecisionTYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
	                PrecisionTYPE *gradientValuesZ = &gradientValuesY[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
	                PrecisionTYPE newGradientValueX, newGradientValueY, newGradientValueZ;
	                for(int i=0; i<controlPointImage->nx*controlPointImage->ny*controlPointImage->nz; i++){

                        newGradientValueX = *gradientValuesX * reorient_NMI_gradient.m[0][0] +
                                            *gradientValuesY * reorient_NMI_gradient.m[0][1] +
                                            *gradientValuesZ * reorient_NMI_gradient.m[0][2];
                        newGradientValueY = *gradientValuesX * reorient_NMI_gradient.m[1][0] +
                                            *gradientValuesY * reorient_NMI_gradient.m[1][1] +
                                            *gradientValuesZ * reorient_NMI_gradient.m[1][2];
                        newGradientValueZ = *gradientValuesX * reorient_NMI_gradient.m[2][0] +
                                            *gradientValuesY * reorient_NMI_gradient.m[2][1] +
                                            *gradientValuesZ * reorient_NMI_gradient.m[2][2];
                        *gradientValuesX++ = newGradientValueX;
                        *gradientValuesY++ = newGradientValueY;
                        *gradientValuesZ++ = newGradientValueZ;
	                }
                }
                if(flag->gradientSmoothingFlag){
                    reg_gaussianSmoothing<PrecisionTYPE>(   nodeNMIGradientImage,
                                                         param->gradientSmoothingValue,
                                                         NULL);
                }

                /* The other gradients are calculated */
                if(flag->beGradFlag && flag->bendingEnergyFlag && param->bendingEnergyWeight>0){
					reg_bspline_bendingEnergyGradient<PrecisionTYPE>(   controlPointImage,
																	    targetImage,
																	    nodeNMIGradientImage,
																	    param->bendingEnergyWeight);
				}
				if(flag->jlGradFlag && flag->jacobianWeightFlag && param->jacobianWeight>0){
                    reg_bspline_jacobianDeterminantGradient<PrecisionTYPE>( controlPointImage,
                                                                            targetImage,
                                                                            nodeNMIGradientImage,
                                                                            param->jacobianWeight,
                                                                            flag->appJacobianFlag);
				}

				/* The conjugate gradient is computed */
				if(!flag->noConjugateGradient){
					if(iteration==1){
						// first conjugate gradient iteration
						if(flag->twoDimRegistration){
							PrecisionTYPE *conjGPtrX = &conjugateG[0];
							PrecisionTYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
							PrecisionTYPE *conjHPtrX = &conjugateH[0];
							PrecisionTYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
							PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
							PrecisionTYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny];
							for(int i=0; i<nodeNMIGradientImage->nx*nodeNMIGradientImage->ny;i++){
								*conjHPtrX++ = *conjGPtrX++ = - *gradientValuesX++;
								*conjHPtrY++ = *conjGPtrY++ = - *gradientValuesY++;
							}
						}else{
							PrecisionTYPE *conjGPtrX = &conjugateG[0];
							PrecisionTYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *conjGPtrZ = &conjGPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *conjHPtrX = &conjugateH[0];
							PrecisionTYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *conjHPtrZ = &conjHPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
							PrecisionTYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
							PrecisionTYPE *gradientValuesZ = &gradientValuesY[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
							for(int i=0; i<nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz;i++){
								*conjHPtrX++ = *conjGPtrX++ = - *gradientValuesX++;
								*conjHPtrY++ = *conjGPtrY++ = - *gradientValuesY++;
								*conjHPtrZ++ = *conjGPtrZ++ = - *gradientValuesZ++;
							}
						}
					}
					else{
						double dgg=0.0, gg=0.0;
						if(flag->twoDimRegistration){
							PrecisionTYPE *conjGPtrX = &conjugateG[0];
							PrecisionTYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
							PrecisionTYPE *conjHPtrX = &conjugateH[0];
							PrecisionTYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny];
							PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
							PrecisionTYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny];
							for(int i=0; i<nodeNMIGradientImage->nx*nodeNMIGradientImage->ny;i++){
								gg += conjHPtrX[i] * conjGPtrX[i];
								gg += conjHPtrY[i] * conjGPtrY[i];
								dgg += (gradientValuesX[i] + conjGPtrX[i]) * gradientValuesX[i];
								dgg += (gradientValuesY[i] + conjGPtrY[i]) * gradientValuesY[i];
							}
							double gam = dgg/gg;
							for(int i=0; i<nodeNMIGradientImage->nx*nodeNMIGradientImage->ny;i++){
								conjGPtrX[i] = - gradientValuesX[i];
								conjGPtrY[i] = - gradientValuesY[i];
								conjHPtrX[i] = (float)(conjGPtrX[i] + gam * conjHPtrX[i]);
								conjHPtrY[i] = (float)(conjGPtrY[i] + gam * conjHPtrY[i]);
								gradientValuesX[i] = - conjHPtrX[i];
								gradientValuesY[i] = - conjHPtrY[i];
							}
						}
						else{
							PrecisionTYPE *conjGPtrX = &conjugateG[0];
							PrecisionTYPE *conjGPtrY = &conjGPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *conjGPtrZ = &conjGPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *conjHPtrX = &conjugateH[0];
							PrecisionTYPE *conjHPtrY = &conjHPtrX[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *conjHPtrZ = &conjHPtrY[nodeNMIGradientImage->nx * nodeNMIGradientImage->ny * nodeNMIGradientImage->nz];
							PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
							PrecisionTYPE *gradientValuesY = &gradientValuesX[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
							PrecisionTYPE *gradientValuesZ = &gradientValuesY[nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz];
							for(int i=0; i<nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz;i++){
								gg += conjHPtrX[i] * conjGPtrX[i];
								gg += conjHPtrY[i] * conjGPtrY[i];
								gg += conjHPtrZ[i] * conjGPtrZ[i];
								dgg += (gradientValuesX[i] + conjGPtrX[i]) * gradientValuesX[i];
								dgg += (gradientValuesY[i] + conjGPtrY[i]) * gradientValuesY[i];
								dgg += (gradientValuesZ[i] + conjGPtrZ[i]) * gradientValuesZ[i];
							}
							double gam = dgg/gg;
							for(int i=0; i<nodeNMIGradientImage->nx*nodeNMIGradientImage->ny*nodeNMIGradientImage->nz;i++){
								conjGPtrX[i] = - gradientValuesX[i];
								conjGPtrY[i] = - gradientValuesY[i];
								conjGPtrZ[i] = - gradientValuesZ[i];
								conjHPtrX[i] = (float)(conjGPtrX[i] + gam * conjHPtrX[i]);
								conjHPtrY[i] = (float)(conjGPtrY[i] + gam * conjHPtrY[i]);
								conjHPtrZ[i] = (float)(conjGPtrZ[i] + gam * conjHPtrZ[i]);
								gradientValuesX[i] = - conjHPtrX[i];
								gradientValuesY[i] = - conjHPtrY[i];
								gradientValuesZ[i] = - conjHPtrZ[i];
							}
						}
					}
				}
				maxLength = reg_getMaximalLength<PrecisionTYPE>(nodeNMIGradientImage);


			/* The gradient is applied to the control point positions */
			if(maxLength==0){
				printf("No Gradient ... exit\n");
				break;	
			}

			/* ** LINE ASCENT ** */
			int lineIteration = 0;
    		currentSize=maxStepSize;
			float addedStep=0.0f;

            while(currentSize>smallestSize && lineIteration<12){

				float currentLength = -currentSize/maxLength;

                  /* Update the control point position */
                    if(flag->useCompositionFlag){ // the control point positions are updated using composition
                        memcpy(controlPointImage->data,bestControlPointPosition,controlPointImage->nvox*controlPointImage->nbyper);
                        reg_spline_cppComposition(  controlPointImage,
                                                    nodeNMIGradientImage,
                                                    (float)currentLength,
                                                    1);
                    }
                    else{ // the control point positions are updated using addition
                        if(flag->twoDimRegistration){
                            PrecisionTYPE *controlPointValuesX = NULL;
                            PrecisionTYPE *controlPointValuesY = NULL;
                            controlPointValuesX = static_cast<PrecisionTYPE *>(controlPointImage->data);
                            controlPointValuesY = &controlPointValuesX[controlPointImage->nx*controlPointImage->ny];
                            PrecisionTYPE *bestControlPointValuesX = &bestControlPointPosition[0];
                            PrecisionTYPE *bestControlPointValuesY = &bestControlPointValuesX[controlPointImage->nx*controlPointImage->ny];
                            PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
                            PrecisionTYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx*controlPointImage->ny];
                            for(int i=0; i<controlPointImage->nx*controlPointImage->ny;i++){
                                *controlPointValuesX++ = *bestControlPointValuesX++ + currentLength * *gradientValuesX++;
                                *controlPointValuesY++ = *bestControlPointValuesY++ + currentLength * *gradientValuesY++;
                            }
                        }
                        else{
                            PrecisionTYPE *controlPointValuesX = NULL;
                            PrecisionTYPE *controlPointValuesY = NULL;
                            PrecisionTYPE *controlPointValuesZ = NULL;
                            controlPointValuesX = static_cast<PrecisionTYPE *>(controlPointImage->data);
                            controlPointValuesY = &controlPointValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                            controlPointValuesZ = &controlPointValuesY[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                            PrecisionTYPE *bestControlPointValuesX = &bestControlPointPosition[0];
                            PrecisionTYPE *bestControlPointValuesY = &bestControlPointValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                            PrecisionTYPE *bestControlPointValuesZ = &bestControlPointValuesY[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                            PrecisionTYPE *gradientValuesX = static_cast<PrecisionTYPE *>(nodeNMIGradientImage->data);
                            PrecisionTYPE *gradientValuesY = &gradientValuesX[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                            PrecisionTYPE *gradientValuesZ = &gradientValuesY[controlPointImage->nx*controlPointImage->ny*controlPointImage->nz];
                            for(int i=0; i<controlPointImage->nx*controlPointImage->ny*controlPointImage->nz;i++){
                                *controlPointValuesX++ = *bestControlPointValuesX++ + currentLength * *gradientValuesX++;
                                *controlPointValuesY++ = *bestControlPointValuesY++ + currentLength * *gradientValuesY++;
                                *controlPointValuesZ++ = *bestControlPointValuesZ++ + currentLength * *gradientValuesZ++;
                            }
                        }
                    }

                /* The Jacobian-based penalty term is computed */
                    if(flag->jacobianWeightFlag && param->jacobianWeight>0){
                        currentWJac = param->jacobianWeight
                            * reg_bspline_jacobian<PrecisionTYPE>(controlPointImage, targetImage, flag->appJacobianFlag);
                    }

                int negCorrection=0;
                while(currentWJac!=currentWJac && negCorrection<5){
                    printf("*");
                        currentWJac = param->jacobianWeight*
                            reg_bspline_correctFolding<PrecisionTYPE>(  controlPointImage,
                                                                        targetImage,
                                                                        flag->appJacobianFlag);
                    negCorrection++;
                }

                /* The bending energy is computed */
                    if(flag->bendingEnergyFlag && param->bendingEnergyWeight>0){
                        currentWBE = param->bendingEnergyWeight
                                * reg_bspline_bendingEnergy<PrecisionTYPE>(controlPointImage, targetImage, flag->appBendingEnergyFlag);
                    }


                /* The source image is resampled and the metric evaluated */
					reg_bspline<PrecisionTYPE>(	controlPointImage,
											    targetImage,
											    positionFieldImage,
											    targetMask,
											    0);

					/* Resample the source image */
					reg_resampleSourceImage<PrecisionTYPE>(	targetImage,
										                    sourceImage,
										                    resultImage,
										                    positionFieldImage,
                                                            NULL,
										                    1,
										                    param->sourcePaddingValue);
				/* Computation of the Metric Value */
				if(flag->useSSDFlag){
					SSDValue=reg_getSSD<PrecisionTYPE>(	targetImage,
                                                        resultImage);
                    currentValue = -(1.0-param->bendingEnergyWeight-param->jacobianWeight)*log(SSDValue+1.0);
				}
				else{
                    reg_getEntropies<double>(	targetImage,
											    resultImage,
											    JH_PW_APPROX,
											    param->binning,
											    probaJointHistogram,
											    logJointHistogram,
											    entropies,
                                                targetMask);
                    currentValue = (1.0-param->bendingEnergyWeight-param->jacobianWeight)*(entropies[0]+entropies[1])/entropies[2];
                }

                currentValue -= currentWBE + currentWJac;

                iteration++;
                lineIteration++;

				/* The current deformation field is kept if it was the best so far */
				if(currentValue>bestValue){
						memcpy(bestControlPointPosition,controlPointImage->data,controlPointImage->nvox*controlPointImage->nbyper);
					bestValue = currentValue;
					bestWBE = currentWBE;
					bestWJac = currentWJac;
					addedStep += currentSize;
					currentSize*=1.1f;
					currentSize = (currentSize<maxStepSize)?currentSize:maxStepSize;
				}
				else{
					currentSize*=0.5;
				}
			}
				memcpy(controlPointImage->data,bestControlPointPosition,controlPointImage->nvox*controlPointImage->nbyper);
			currentSize=addedStep;
			printf("[%i] Objective function value=%g | max. added disp. = %g mm", iteration, bestValue, addedStep);
			if(flag->bendingEnergyFlag && param->bendingEnergyWeight>0) printf(" | wBE=%g", bestWBE);
			if(flag->jacobianWeightFlag && param->jacobianWeight>0) printf(" | wJacLog=%g", bestWJac);
			printf("\n");
		} // while(iteration<param->maxIteration && currentSize>smallestSize){

        /* The deformation model is unfolded if the Jacobian
            determinant-based penalty term is used */
        if(flag->jacobianWeightFlag && param->jacobianWeight>0){
            int finalNegCorrection=0;
            PrecisionTYPE finalWJac=0;
            do{
                {
                    finalWJac = param->jacobianWeight*
                        reg_bspline_correctFolding<PrecisionTYPE>(controlPointImage,
                                                                  targetImage,
                                                                  false); // No approximation is done
                }
                finalNegCorrection++;
                if(finalWJac!=finalWJac)
                    printf( "*** Final folding correction [%i/%i] ***\n",
                        finalNegCorrection, FOLDING_CORRECTION_STEP);
                else printf(">>> Final Jacobian based penalty term value = %g\n", finalWJac);
            }
            while(finalWJac!=finalWJac && finalNegCorrection<FOLDING_CORRECTION_STEP);

            if(finalWJac!=finalWJac){
                fprintf(stderr, "[WARNING] The final folding correction failed.\n");
                fprintf(stderr, "[WARNING] You might want to increase the penalty term weights.\n");
                {
                    memcpy(controlPointImage->data,bestControlPointPosition,controlPointImage->nvox*controlPointImage->nbyper);
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
			if(!flag->noConjugateGradient){
				free(conjugateG);
				free(conjugateH);
            }
            free(bestControlPointPosition);

        nifti_image_free( resultImage );


        /* ****************** */
        /* OUTPUT THE RESULTS */
        /* ****************** */
        if(level==(param->level2Perform-1)){

                if(param->level2Perform != param->levelNumber){
                    if(positionFieldImage->data)free(positionFieldImage->data);
                    positionFieldImage->dim[1]=positionFieldImage->nx=targetHeader->nx;
                    positionFieldImage->dim[2]=positionFieldImage->ny=targetHeader->ny;
                    positionFieldImage->dim[3]=positionFieldImage->nz=targetHeader->nz;
                    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
                    if(flag->twoDimRegistration)
                        positionFieldImage->dim[5]=positionFieldImage->nu=2;
                    else positionFieldImage->dim[5]=positionFieldImage->nu=3;
                    positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
                    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
                    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
                    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
                    positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
                }

			/* The corresponding deformation field is evaluated and saved */
			/* The best result is returned */
            nifti_set_filenames(controlPointImage, param->outputCPPName, 0, 0);
			nifti_image_write(controlPointImage);

			reg_bspline<PrecisionTYPE>(	controlPointImage,
										targetHeader,
										positionFieldImage,
										NULL,
										0);

            nifti_image_free( sourceImage );
            sourceImage = nifti_image_read(param->sourceImageName,true); // reload the source image with the correct intensity values

			resultImage = nifti_copy_nim_info(targetHeader);
			resultImage->cal_min = sourceImage->cal_min;
			resultImage->cal_max = sourceImage->cal_max;
			resultImage->scl_slope = sourceImage->scl_slope;
			resultImage->scl_inter = sourceImage->scl_inter;
			resultImage->datatype = sourceImage->datatype;
			resultImage->nbyper = sourceImage->nbyper;
			resultImage->data = (void *)calloc(resultImage->nvox, resultImage->nbyper);
			reg_resampleSourceImage<double>(targetHeader,
                                            sourceImage,
							                resultImage,
							                positionFieldImage,
                                            NULL,
							                3,
							                param->sourcePaddingValue);
			if(!flag->outputResultFlag) param->outputResultName=(char *)"outputResult.nii";
			nifti_set_filenames(resultImage, param->outputResultName, 0, 0);
			nifti_image_write(resultImage);
			nifti_image_free(resultImage);

		} // if(level==(param->levelNumber-1)){

		nifti_image_free( positionFieldImage );
		nifti_image_free( sourceImage );
		nifti_image_free( targetImage );

		printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
	} // for(int level=0; level<param->levelNumber; level++){

	/* Mr Clean */
	nifti_image_free( controlPointImage );
	nifti_image_free( targetHeader );
	nifti_image_free( sourceHeader );
    if(flag->targetMaskFlag)
		nifti_image_free( targetMaskImage );

	free( flag );
	free( param );

	time_t end; time( &end );
	int minutes = (int)floorf(float(end-start)/60.0f);
	int seconds = (int)(end-start - 60*minutes);
	printf("Registration Performed in %i min %i sec\n", minutes, seconds);
	printf("Have a good day !\n");

	return 0;
}
