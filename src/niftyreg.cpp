// Registration types (degrees of freedom) - these constants have to have these values
#define RIGID 0
#define AFFINE 1

// Precision levels for final interpolation
#define INTERP_PREC_SOURCE 0
#define INTERP_PREC_DOUBLE 1

// Working precision for the source and target images: float and double are valid for the NiftyReg code
#define PRECISION_TYPE double

#include "_reg_resampling.h"
#include "_reg_tools.h"
#include "_reg_aladin.h"
#include "_reg_f3d.h"
#include "_reg_f3d_sym.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "niftyreg.h"
#include "substitutions.h"

extern "C"
SEXP reg_aladin_R (SEXP source, SEXP target, SEXP type, SEXP finalPrecision, SEXP nLevels, SEXP maxIterations, SEXP useBlockPercentage, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP verbose)
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
    free(result.completedIterations);
    
    UNPROTECT(affineProvided ? 3 : 4);
    
    return returnValue;
}

extern "C"
SEXP reg_f3d_R (SEXP source, SEXP target, SEXP finalPrecision, SEXP nLevels, SEXP maxIterations, SEXP nBins, SEXP bendingEnergyWeight, SEXP jacobianWeight, SEXP finalSpacing, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP initControl, SEXP verbose)
{
    int i, j, levels = *(INTEGER(nLevels));
    SEXP returnValue, data, controlPoints, completedIterations, xform;
    
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
    
    float spacing[3];
    for (i = 0; i < 3; i++)
        spacing[i] = (float) REAL(finalSpacing)[i];
    
    int precisionType = (strcmp(CHAR(STRING_ELT(finalPrecision,0)),"source")==0 ? INTERP_PREC_SOURCE : INTERP_PREC_DOUBLE);
    
    f3d_result result = do_reg_f3d(sourceImage, targetImage, precisionType, *(INTEGER(nLevels)), *(INTEGER(maxIterations)), *(INTEGER(finalInterpolation)), NULL, targetMaskImage, controlPointImage, affineTransformation, *(INTEGER(nBins)), spacing, (float) *(REAL(bendingEnergyWeight)), (float) *(REAL(jacobianWeight)), 0.01, false, (*(INTEGER(verbose)) == 1));
    
    PROTECT(returnValue = NEW_LIST(4));
    
    if (result.forwardImage->datatype == DT_INT32)
    {
        // Integer-valued data went in, and precision must be "source"
        PROTECT(data = NEW_INTEGER((R_len_t) result.forwardImage->nvox));
        for (size_t i = 0; i < result.forwardImage->nvox; i++)
            INTEGER(data)[i] = ((int *) result.forwardImage->data)[i];
    }
    else
    {
        // All other cases
        PROTECT(data = NEW_NUMERIC((R_len_t) result.forwardImage->nvox));
        for (size_t i = 0; i < result.forwardImage->nvox; i++)
            REAL(data)[i] = ((double *) result.forwardImage->data)[i];
    }
    
    SET_ELEMENT(returnValue, 0, data);
    
    PROTECT(controlPoints = NEW_NUMERIC((R_len_t) result.forwardControlPoints->nvox));
    for (size_t i = 0; i < result.forwardControlPoints->nvox; i++)
        REAL(controlPoints)[i] = ((double *) result.forwardControlPoints->data)[i];
    
    SET_ELEMENT(returnValue, 1, controlPoints);
    
    PROTECT(completedIterations = NEW_INTEGER((R_len_t) levels));
    for (i = 0; i < levels; i++)
        INTEGER(completedIterations)[i] = result.completedIterations[i];
    
    SET_ELEMENT(returnValue, 2, completedIterations);
    
    PROTECT(xform = NEW_NUMERIC(21));
    REAL(xform)[0] = (double) result.forwardControlPoints->qform_code;
    REAL(xform)[1] = (double) result.forwardControlPoints->sform_code;
    REAL(xform)[2] = (double) result.forwardControlPoints->quatern_b;
    REAL(xform)[3] = (double) result.forwardControlPoints->quatern_c;
    REAL(xform)[4] = (double) result.forwardControlPoints->quatern_d;
    REAL(xform)[5] = (double) result.forwardControlPoints->qoffset_x;
    REAL(xform)[6] = (double) result.forwardControlPoints->qoffset_y;
    REAL(xform)[7] = (double) result.forwardControlPoints->qoffset_z;
    REAL(xform)[8] = (double) result.forwardControlPoints->qfac;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
            REAL(xform)[(i*4)+j+9] = (double) result.forwardControlPoints->sto_xyz.m[i][j];
    }
    
    SET_ELEMENT(returnValue, 3, xform);
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    if (targetMaskImage != NULL)
        nifti_image_free(targetMaskImage);
    if (controlPointImage != NULL)
        nifti_image_free(controlPointImage);
    nifti_image_free(result.forwardImage);
    nifti_image_free(result.forwardControlPoints);
    free(result.completedIterations);
    if (affineProvided || isNull(initControl))
        free(affineTransformation);
    
    UNPROTECT(5);
    
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

mat44 * create_init_affine (nifti_image *sourceImage, nifti_image *targetImage)
{
    mat44 *affineTransformation = (mat44 *) calloc(1, sizeof(mat44));
    for (int i = 0; i < 4; i++)
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
	
    return affineTransformation;
}

// Run the "aladin" registration algorithm
aladin_result do_reg_aladin (nifti_image *sourceImage, nifti_image *targetImage, int type, int finalPrecision, int nLevels, int maxIterations, int useBlockPercentage, int finalInterpolation, nifti_image *targetMaskImage, mat44 *affineTransformation, bool verbose)
{
    if (affineTransformation == NULL)
        affineTransformation = create_init_affine(sourceImage, targetImage);
    
    reg_aladin<PRECISION_TYPE> *reg = new reg_aladin<PRECISION_TYPE>;
    
    reg->SetMaxIterations(maxIterations);
    reg->SetNumberOfLevels(nLevels);
    reg->SetLevelsToPerform(nLevels);
    reg->SetReferenceSigma(0.0);
    reg->SetFloatingSigma(0.0);
    reg->SetAlignCentre(1);
    reg->SetPerformAffine(type == AFFINE);
    reg->SetPerformRigid(1);
    reg->SetVerbose((int) verbose);
    reg->SetBlockPercentage(useBlockPercentage);
    reg->SetInlierLts(50.0);
    reg->SetInterpolation(finalInterpolation);
    
    // Set the reference and floating images
    reg->SetInputReference(targetImage);
    reg->SetInputFloating(sourceImage);
    
    // Set the initial affine transformation
    reg->SetTransformationMatrix(affineTransformation);
    
    // Set the reference mask if defined
    if (targetMaskImage != NULL)
        reg->SetInputMask(targetMaskImage);
    
    // Run the registration
    reg->Run();
    
    int *completedIterations = (int *) calloc(nLevels, sizeof(int));
    memcpy(completedIterations, reg->GetCompletedIterations(), nLevels*sizeof(int));
    
    // Store the results
    aladin_result result;
    result.image = copy_complete_nifti_image(reg->GetFinalWarpedImage());
    result.affine = (mat44 *) calloc(1, sizeof(mat44));
    memcpy(result.affine, reg->GetTransformationMatrix(), sizeof(mat44));
    result.completedIterations = completedIterations;
    
    delete reg;
    
    return result;
}

f3d_result do_reg_f3d (nifti_image *sourceImage, nifti_image *targetImage, int finalPrecision, int nLevels, int maxIterations, int finalInterpolation, nifti_image *sourceMaskImage, nifti_image *targetMaskImage, nifti_image *controlPointImage, mat44 *affineTransformation, int nBins, float *spacing, float bendingEnergyWeight, float jacobianWeight, float inverseConsistencyWeight, bool symmetric, bool verbose)
{
    if (controlPointImage == NULL && affineTransformation == NULL)
        affineTransformation = create_init_affine(sourceImage, targetImage);
    
    // Binarise the mask image
    if (targetMaskImage != NULL)
        reg_tools_binarise_image(targetMaskImage);
    
    // The source data type is changed for precision if requested
    if (finalPrecision == INTERP_PREC_DOUBLE)
        reg_tools_changeDatatype<double>(sourceImage);
    
    f3d_result result;
    
    if (nLevels == 0)
    {
        // Allocate and set completed iterations (nothing will be done so this is zero)
        int *completedIterations = (int *) calloc(1, sizeof(int));
        *completedIterations = 0;
        
        reg_checkAndCorrectDimension(controlPointImage);
        
        // Allocate a deformation field image
        nifti_image *deformationFieldImage = create_position_field(targetImage, targetImage->nz <= 1);
        
        // Calculate deformation field from the control point image
        if(controlPointImage->pixdim[5] > 1)
            reg_bspline_getDeformationFieldFromVelocityGrid(controlPointImage, deformationFieldImage);
        else
            reg_spline_getDeformationField(controlPointImage, targetImage, deformationFieldImage, NULL, false, true);
        
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
        resultImage->nvox = resultImage->dim[1] * resultImage->dim[2] * resultImage->dim[3] * resultImage->dim[4];
        resultImage->data = (void *) calloc(resultImage->nvox, resultImage->nbyper);
        
        // Resample source image to target space
        reg_resampleSourceImage(targetImage, sourceImage, resultImage, deformationFieldImage, NULL, finalInterpolation, 0);
        
        // Free deformation field image
        nifti_image_free(deformationFieldImage);
        
        result.forwardImage = resultImage;
        result.forwardControlPoints = copy_complete_nifti_image(controlPointImage);
        result.reverseImage = NULL;
        result.reverseControlPoints = NULL;
        result.completedIterations = completedIterations;
    }
    else
    {
        reg_f3d<PRECISION_TYPE> *reg = NULL;

        // Create the reg_f3d object
        if (symmetric)
            reg = new reg_f3d_sym<PRECISION_TYPE>(targetImage->nt, sourceImage->nt);
        else
            reg = new reg_f3d<PRECISION_TYPE>(targetImage->nt, sourceImage->nt);
        
#ifdef _OPENMP
        int maxThreadNumber = omp_get_max_threads();
        if (verbose)
            Rprintf("[NiftyReg F3D] Using OpenMP with %i thread(s)\n", maxThreadNumber);
#endif

        // Set the reg_f3d parameters
        reg->SetReferenceImage(targetImage);
        reg->SetFloatingImage(sourceImage);
        
        if (verbose)
            reg->PrintOutInformation();
        else
            reg->DoNotPrintOutInformation();
        
        if (targetMaskImage != NULL)
           reg->SetReferenceMask(targetMaskImage);
        
        if (controlPointImage != NULL)
           reg->SetControlPointGridImage(controlPointImage);
        
        if (affineTransformation != NULL)
           reg->SetAffineTransformation(affineTransformation);
        
        reg->SetBendingEnergyWeight(bendingEnergyWeight);
        reg->SetLinearEnergyWeights(0.0, 0.0);
        reg->SetL2NormDisplacementWeight(0.0);
        reg->SetJacobianLogWeight(jacobianWeight);
        
        reg->SetMaximalIterationNumber(maxIterations);
        
        for (int i = 0; i < 3; i++)
            reg->SetSpacing((unsigned) i, (PRECISION_TYPE) spacing[i]);
        
        reg->SetLevelNumber(nLevels);
        reg->SetLevelToPerform(nLevels);

        if (finalInterpolation == 3)
            reg->UseCubicSplineInterpolation();
        else if (finalInterpolation == 1)
            reg->UseLinearInterpolation();
        else
            reg->UseNeareatNeighborInterpolation();
        
        // Parameters only relevant to the symmetric version of the algorithm
        if (symmetric)
        {
            if (sourceMaskImage != NULL)
                reg->SetFloatingMask(sourceMaskImage);
            
            reg->SetInverseConsistencyWeight(inverseConsistencyWeight);
        }
        
        // Run the registration
        reg->Run_f3d();
        
        int *completedIterations = (int *) calloc(nLevels, sizeof(int));
        memcpy(completedIterations, reg->GetCompletedIterations(), nLevels*sizeof(int));
        
        result.forwardImage = reg->GetWarpedImage()[0];
        result.forwardControlPoints = reg->GetControlPointPositionImage();
        if (symmetric)
        {
            result.reverseImage = reg->GetWarpedImage()[1];
            result.reverseControlPoints = reg->GetBackwardControlPointPositionImage();
        }
        result.completedIterations = completedIterations;
        
        // Erase the registration object
        delete reg;
    }
    
    return result;
}
