#include <RcppEigen.h>

#include "RNifti.h"

#include "helpers.h"
#include "_reg_tools.h"

using namespace RNifti;

int nonunitaryDims (const NiftiImage &image)
{
    int nDims = image.nDims();
    while (nDims > 0 && image->dim[nDims] < 2)
        nDims--;
    return nDims;
}

bool isMultichannel (const NiftiImage &image)
{
    // Assume 2D RGB or RGBA image
    return (image.nDims() == 3 && (image->nz == 3 || image->nz == 4));
}

NiftiImage collapseChannels (const NiftiImage &image)
{
    if (isMultichannel(image))
    {
        std::vector<double> red = image.slice(0).getData<double>();
        const std::vector<double> green = image.slice(1).getData<double>();
        const std::vector<double> blue = image.slice(2).getData<double>();
        
        for (size_t i=0; i<red.size(); i++)
            red[i] = (red[i] + green[i] + blue[i]) / 3.0;
        
        nifti_image *result = nifti_copy_nim_info(image);
        result->dim[0] = image->dim[0] - 1;
        result->dim[image->dim[0]] = 1;
        result->pixdim[image->dim[0]] = 1.0;
        nifti_update_dims_from_array(result);
        
        result->datatype = DT_FLOAT64;
        nifti_datatype_sizes(result->datatype, &result->nbyper, &result->swapsize);
        
        result->data = calloc(result->nvox, 8);
        std::copy(red.begin(), red.end(), static_cast<double *>(result->data));
        
        return NiftiImage(result);
    }
    else
        return image;
}

void checkImages (const NiftiImage &sourceImage, const NiftiImage &targetImage)
{
    if (sourceImage.isNull())
        throw std::runtime_error("Cannot read or retrieve source image");
    if (targetImage.isNull())
        throw std::runtime_error("Cannot read or retrieve target image");
    
    // Find the dimensionality of the images, ignoring final unitary dimensions
    const int nSourceDim = nonunitaryDims(sourceImage);
    const int nTargetDim = nonunitaryDims(targetImage);
    
    if (nSourceDim < 2 || nSourceDim > 4)
        throw std::runtime_error("Source image should have 2, 3 or 4 dimensions");
    if (nTargetDim < 2 || nTargetDim > 3)
        throw std::runtime_error("Target image should have 2 or 3 dimensions");
    
    const std::vector<int> sourceDims = sourceImage.dim();
    const std::vector<int> targetDims = targetImage.dim();
    
    for (int i=0; i<std::min(nSourceDim,nTargetDim); i++)
    {
        if (sourceDims[i] < 4 && (i < (nSourceDim-1) || !isMultichannel(sourceImage)))
            throw std::runtime_error("Source image should have width at least 4 in all dimensions");
    }
    for (int i=0; i<nTargetDim; i++)
    {
        if (targetDims[i] < 4 && (i < (nTargetDim-1) || !isMultichannel(targetImage)))
            throw std::runtime_error("Target image should have width at least 4 in all dimensions");
    }
}

NiftiImage normaliseImage (const NiftiImage &image)
{
    NiftiImage normalisedImage = image;
    normalisedImage.drop();
    
    // NiftyReg assumes RAS when qform and sform codes are both zero
    // That's reasonable for NIfTI, but ANALYZE defaults to LAS
    if (image->qform_code == 0 && image->sform_code == 0 && image->nifti_type == NIFTI_FTYPE_ANALYZE)
    {
        Rf_warning("ANALYZE voxel order will be assumed to be LAS, but use of this format is strongly discouraged");
        mat44 defaultXform = image.xform();
        defaultXform.m[0][0] *= -1.0;
        normalisedImage->qto_xyz = defaultXform;
        normalisedImage->qto_ijk = nifti_mat44_inverse(defaultXform);
        nifti_mat44_to_quatern(normalisedImage->qto_xyz, &normalisedImage->quatern_b, &normalisedImage->quatern_c, &normalisedImage->quatern_d, &normalisedImage->qoffset_x, &normalisedImage->qoffset_y, &normalisedImage->qoffset_z, NULL, NULL, NULL, &normalisedImage->qfac);
        normalisedImage->qform_code = 2;
    }
    
    // Get everything else the way NiftyReg wants it
    reg_checkAndCorrectDimension(normalisedImage);
    
    return normalisedImage;
}

NiftiImage allocateMultiregResult (const NiftiImage &source, const NiftiImage &target, const bool forceDouble)
{
    nifti_image *newStruct = nifti_copy_nim_info(target);
    newStruct->dim[0] = source->dim[0];
    newStruct->dim[source.nDims()] = source->dim[source.nDims()];
    newStruct->pixdim[source.nDims()] = source->pixdim[source.nDims()];
    
    if (forceDouble)
    {
        newStruct->datatype = DT_FLOAT64;
        nifti_datatype_sizes(newStruct->datatype, &newStruct->nbyper, NULL);
    }
    
    nifti_update_dims_from_array(newStruct);
    
    size_t dataSize = nifti_get_volsize(newStruct);
    newStruct->data = calloc(1, dataSize);
    
    return NiftiImage(newStruct);
}
