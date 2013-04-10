// Registration types (degrees of freedom)
#define TYPE_RIGID 0
#define TYPE_AFFINE 1

// Working precision: float and double are valid for the NiftyReg code
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
SEXP reg_aladin_R (SEXP source, SEXP target, SEXP type, SEXP nLevels, SEXP maxIterations, SEXP useBlockPercentage, SEXP finalInterpolation, SEXP targetMask, SEXP affineComponents, SEXP verbose, SEXP estimateOnly)
{
    int i, j, levels = *(INTEGER(nLevels));
    SEXP returnValue, completedIterations;
    
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
    
    int regType = (strcmp(CHAR(STRING_ELT(type,0)),"rigid")==0 ? TYPE_RIGID : TYPE_AFFINE);
    
    aladin_result result = do_reg_aladin(sourceImage, targetImage, regType, *(INTEGER(nLevels)), *(INTEGER(maxIterations)), *(INTEGER(useBlockPercentage)), *(INTEGER(finalInterpolation)), targetMaskImage, affineTransformation, (*(INTEGER(verbose)) == 1), (*(INTEGER(estimateOnly)) == 1));
    
    PROTECT(returnValue = NEW_LIST(3));
    
    convert_and_insert_image(result.image, returnValue, 0);
    convert_and_insert_affine(result.affine, returnValue, 1);
    
    if (result.completedIterations != NULL)
    {
        PROTECT(completedIterations = NEW_INTEGER((R_len_t) levels));
        for (i = 0; i < levels; i++)
            INTEGER(completedIterations)[i] = result.completedIterations[i];
    
        SET_ELEMENT(returnValue, 2, completedIterations);
        UNPROTECT(1);
    }
    else
        SET_ELEMENT(returnValue, 2, R_NilValue);
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    if (targetMaskImage != NULL)
        nifti_image_free(targetMaskImage);
    if (result.image != NULL)
        nifti_image_free(result.image);
    free(result.affine);
    if (result.completedIterations != NULL)
        free(result.completedIterations);
    
    UNPROTECT(1);
    
    return returnValue;
}

extern "C"
SEXP reg_f3d_R (SEXP source, SEXP target, SEXP nLevels, SEXP maxIterations, SEXP nBins, SEXP bendingEnergyWeight, SEXP jacobianWeight, SEXP inverseConsistencyWeight, SEXP finalSpacing, SEXP finalInterpolation, SEXP targetMask, SEXP sourceMask, SEXP affineComponents, SEXP initControl, SEXP symmetric, SEXP verbose, SEXP estimateOnly)
{
    int i, j, levels = *(INTEGER(nLevels));
    SEXP returnValue, completedIterations;
    
    bool affineProvided = !isNull(affineComponents);
    bool useSymmetricAlgorithm = (*(INTEGER(symmetric)) == 1);
    
    nifti_image *sourceImage = s4_image_to_struct(source);
    nifti_image *targetImage = s4_image_to_struct(target);
    nifti_image *sourceMaskImage = NULL;
    nifti_image *targetMaskImage = NULL;
    nifti_image *controlPointImage = NULL;
    mat44 *affineTransformation = NULL;
    
    if (!isNull(sourceMask) && useSymmetricAlgorithm)
        sourceMaskImage = s4_image_to_struct(sourceMask);
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
        
    f3d_result result = do_reg_f3d(sourceImage, targetImage, *(INTEGER(nLevels)), *(INTEGER(maxIterations)), *(INTEGER(finalInterpolation)), sourceMaskImage, targetMaskImage, controlPointImage, affineTransformation, *(INTEGER(nBins)), spacing, (float) *(REAL(bendingEnergyWeight)), (float) *(REAL(jacobianWeight)), (float) *(REAL(inverseConsistencyWeight)), useSymmetricAlgorithm, (*(INTEGER(verbose)) == 1), (*(INTEGER(estimateOnly)) == 1));
    
    if (useSymmetricAlgorithm)
        PROTECT(returnValue = NEW_LIST(8));
    else
        PROTECT(returnValue = NEW_LIST(4));
    
    convert_and_insert_image(result.forwardImage, returnValue, 0);
    convert_and_insert_image(result.forwardControlPoints, returnValue, 1);
    convert_and_insert_xform(result.forwardControlPoints, returnValue, 2);
    
    if (result.completedIterations != NULL)
    {
        PROTECT(completedIterations = NEW_INTEGER((R_len_t) levels));
        for (i = 0; i < levels; i++)
            INTEGER(completedIterations)[i] = result.completedIterations[i];
    
        SET_ELEMENT(returnValue, 3, completedIterations);
        UNPROTECT(1);
    }
    else
        SET_ELEMENT(returnValue, 3, R_NilValue);
    
    if (useSymmetricAlgorithm)
    {
        convert_and_insert_image(result.reverseImage, returnValue, 4);
        convert_and_insert_image(result.reverseControlPoints, returnValue, 5);
        convert_and_insert_xform(result.reverseControlPoints, returnValue, 6);
        convert_and_insert_affine(result.initAffine, returnValue, 7);
    }
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    if (targetMaskImage != NULL)
        nifti_image_free(targetMaskImage);
    if (controlPointImage != NULL)
        nifti_image_free(controlPointImage);
    if (result.forwardImage != NULL)
        nifti_image_free(result.forwardImage);
    nifti_image_free(result.forwardControlPoints);
    if (result.completedIterations != NULL)
        free(result.completedIterations);
    if (affineProvided || isNull(initControl))
        free(affineTransformation);
    
    UNPROTECT(1);
    
    return returnValue;
}

void convert_and_insert_image (nifti_image *image, SEXP list, int index)
{
    if (image == NULL)
    {
        SET_ELEMENT(list, index, R_NilValue);
        return;
    }
    
    SEXP data;
    
    if (image->datatype == DT_INT32)
    {
        PROTECT(data = NEW_INTEGER((R_len_t) image->nvox));
        for (size_t i = 0; i < image->nvox; i++)
            INTEGER(data)[i] = ((int *) image->data)[i];
    }
    else
    {
        PROTECT(data = NEW_NUMERIC((R_len_t) image->nvox));
        for (size_t i = 0; i < image->nvox; i++)
            REAL(data)[i] = ((double *) image->data)[i];
    }
    
    SET_ELEMENT(list, index, data);
    UNPROTECT(1);
}

void convert_and_insert_xform (nifti_image *image, SEXP list, int index)
{
    int i, j;
    SEXP xform;
    
    PROTECT(xform = NEW_NUMERIC(21));
    REAL(xform)[0] = (double) image->qform_code;
    REAL(xform)[1] = (double) image->sform_code;
    REAL(xform)[2] = (double) image->quatern_b;
    REAL(xform)[3] = (double) image->quatern_c;
    REAL(xform)[4] = (double) image->quatern_d;
    REAL(xform)[5] = (double) image->qoffset_x;
    REAL(xform)[6] = (double) image->qoffset_y;
    REAL(xform)[7] = (double) image->qoffset_z;
    REAL(xform)[8] = (double) image->qfac;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 4; j++)
            REAL(xform)[(i*4)+j+9] = (double) image->sto_xyz.m[i][j];
    }
    
    SET_ELEMENT(list, index, xform);
    UNPROTECT(1);
}

void convert_and_insert_affine (mat44 *affine, SEXP list, int index)
{
    int i, j;
    SEXP elements;
    
    PROTECT(elements = NEW_NUMERIC(16));
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
            REAL(elements)[(j*4)+i] = (double) affine->m[i][j];
    }
    
    SET_ELEMENT(list, index, elements);
    UNPROTECT(1);
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

nifti_image * get_deformation_field (nifti_image *targetImage, nifti_image *controlPointImage, mat44 *affineTransformation)
{
    // Allocate a deformation field image
    nifti_image *deformationFieldImage = nifti_copy_nim_info(targetImage);
    
    // Set up properties of deformation field
    deformationFieldImage->dim[0] = deformationFieldImage->ndim = 5;
    deformationFieldImage->dim[1] = deformationFieldImage->nx = targetImage->nx;
    deformationFieldImage->dim[2] = deformationFieldImage->ny = targetImage->ny;
    deformationFieldImage->dim[3] = deformationFieldImage->nz = targetImage->nz;
    deformationFieldImage->dim[4] = deformationFieldImage->nt = 1;
    deformationFieldImage->pixdim[4] = deformationFieldImage->dt = 1.0;
    deformationFieldImage->dim[5] = deformationFieldImage->nu = (targetImage->nz <= 1 ? 2 : 3);
    deformationFieldImage->pixdim[5] = deformationFieldImage->du = 1.0;
    deformationFieldImage->dim[6] = deformationFieldImage->nv = 1;
    deformationFieldImage->pixdim[6] = deformationFieldImage->dv = 1.0;
    deformationFieldImage->dim[7] = deformationFieldImage->nw = 1;
    deformationFieldImage->pixdim[7] = deformationFieldImage->dw = 1.0;
    deformationFieldImage->nvox = deformationFieldImage->nx * deformationFieldImage->ny * deformationFieldImage->nz * deformationFieldImage->nt * deformationFieldImage->nu;
    deformationFieldImage->datatype = (sizeof(PRECISION_TYPE)==4 ? NIFTI_TYPE_FLOAT32 : NIFTI_TYPE_FLOAT64);
    deformationFieldImage->nbyper = sizeof(PRECISION_TYPE);
    deformationFieldImage->data = calloc(deformationFieldImage->nvox, deformationFieldImage->nbyper);
    
    if (controlPointImage != NULL)
    {
        reg_checkAndCorrectDimension(controlPointImage);
        
        // Calculate deformation field from the control point image
        if(controlPointImage->pixdim[5] > 1)
            reg_bspline_getDeformationFieldFromVelocityGrid(controlPointImage, deformationFieldImage);
        else
            reg_spline_getDeformationField(controlPointImage, targetImage, deformationFieldImage, NULL, false, true);
    }
    else
    {
        // Calculate deformation field from the affine matrix
        reg_affine_positionField(affineTransformation, targetImage, deformationFieldImage);
    }
    
    return deformationFieldImage;
}

nifti_image * resample_image (nifti_image *sourceImage, nifti_image *targetImage, nifti_image *controlPointImage, mat44 *affineTransformation, int finalInterpolation)
{
    // Get the deformation field
    nifti_image *deformationFieldImage = get_deformation_field(targetImage, controlPointImage, affineTransformation);
    
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
    
    return resultImage;
}

// Run the "aladin" registration algorithm
aladin_result do_reg_aladin (nifti_image *sourceImage, nifti_image *targetImage, int type, int nLevels, int maxIterations, int useBlockPercentage, int finalInterpolation, nifti_image *targetMaskImage, mat44 *affineTransformation, bool verbose, bool estimateOnly)
{
    if (affineTransformation == NULL)
        affineTransformation = create_init_affine(sourceImage, targetImage);
    
    // Binarise the mask images
    if (targetMaskImage != NULL)
        reg_tools_binarise_image(targetMaskImage);
    
    // The source data type is changed for interpolation precision if necessary
    if (finalInterpolation != 0)
        reg_tools_changeDatatype<double>(sourceImage);
    
    aladin_result result;
    
    if (nLevels == 0)
    {
        result.image = resample_image(sourceImage, targetImage, NULL, affineTransformation, finalInterpolation);
        result.affine = (mat44 *) calloc(1, sizeof(mat44));
        memcpy(result.affine, affineTransformation, sizeof(mat44));
        result.completedIterations = NULL;
    }
    else
    {
        reg_aladin<PRECISION_TYPE> *reg = new reg_aladin<PRECISION_TYPE>;
    
        reg->SetMaxIterations(maxIterations);
        reg->SetNumberOfLevels(nLevels);
        reg->SetLevelsToPerform(nLevels);
        reg->SetReferenceSigma(0.0);
        reg->SetFloatingSigma(0.0);
        reg->SetAlignCentre(1);
        reg->SetPerformAffine((int) type == TYPE_AFFINE);
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
        if (estimateOnly)
            result.image = NULL;
        else
            result.image = copy_complete_nifti_image(reg->GetFinalWarpedImage());
        result.affine = (mat44 *) calloc(1, sizeof(mat44));
        memcpy(result.affine, reg->GetTransformationMatrix(), sizeof(mat44));
        result.completedIterations = completedIterations;
    
        delete reg;
    }
    
    return result;
}

f3d_result do_reg_f3d (nifti_image *sourceImage, nifti_image *targetImage, int nLevels, int maxIterations, int finalInterpolation, nifti_image *sourceMaskImage, nifti_image *targetMaskImage, nifti_image *controlPointImage, mat44 *affineTransformation, int nBins, float *spacing, float bendingEnergyWeight, float jacobianWeight, float inverseConsistencyWeight, bool symmetric, bool verbose, bool estimateOnly)
{
    if (controlPointImage == NULL && affineTransformation == NULL)
        affineTransformation = create_init_affine(sourceImage, targetImage);
    else if (controlPointImage != NULL)
        affineTransformation = NULL;
    
    // Binarise the mask images
    if (targetMaskImage != NULL)
        reg_tools_binarise_image(targetMaskImage);
    if (sourceMaskImage != NULL)
        reg_tools_binarise_image(sourceMaskImage);
    
    // Change data types for interpolation precision if necessary
    if (finalInterpolation != 0)
    {
        reg_tools_changeDatatype<double>(sourceImage);
        if (symmetric)
            reg_tools_changeDatatype<double>(targetImage);
    }
    
    f3d_result result;
    
    if (nLevels == 0)
    {
        result.initAffine = NULL;
        result.forwardImage = resample_image(sourceImage, targetImage, controlPointImage, affineTransformation, finalInterpolation);
        result.forwardControlPoints = copy_complete_nifti_image(controlPointImage);
        result.reverseImage = NULL;
        result.reverseControlPoints = NULL;
        result.completedIterations = NULL;
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
        {
            if (symmetric)
            {
                mat44 inverseTransformation = nifti_mat44_inverse(*affineTransformation);
                if (sourceImage->sform_code > 0)
                    sourceImage->sto_xyz = reg_mat44_mul(&inverseTransformation, &(sourceImage->sto_xyz));
                else
                {
                    sourceImage->sform_code = 1;
                    sourceImage->sto_xyz = reg_mat44_mul(&inverseTransformation, &(sourceImage->qto_xyz));
                }
                sourceImage->sto_ijk = nifti_mat44_inverse(sourceImage->sto_xyz);
            }
            else
                reg->SetAffineTransformation(affineTransformation);
        }
        
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
        
        result.forwardImage = NULL;
        result.reverseImage = NULL;
        
        if (!estimateOnly)
            result.forwardImage = reg->GetWarpedImage()[0];
        result.forwardControlPoints = reg->GetControlPointPositionImage();
        if (symmetric)
        {
            result.initAffine = affineTransformation;
            if (!estimateOnly)
                result.reverseImage = reg->GetWarpedImage()[1];
            result.reverseControlPoints = reg->GetBackwardControlPointPositionImage();
        }
        else
            result.initAffine = NULL;
        result.completedIterations = completedIterations;
        
        // Erase the registration object
        delete reg;
    }
    
    return result;
}
