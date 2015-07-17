#include <RcppEigen.h>

#include "config.h"
#include "NiftiImage.h"
#include "DeformationField.h"
#include "aladin.h"
#include "f3d.h"

// Registration types (degrees of freedom)
#define TYPE_RIGID  0
#define TYPE_AFFINE 1

using namespace Rcpp;

typedef std::vector<float> float_vector;

RcppExport SEXP readNifti (SEXP _file, SEXP _internal)
{
BEGIN_RCPP
    NiftiImage image = retrieveImage(_file);
    
    if (as<bool>(_internal))
        return imageToPointer(image, "NIfTI image");
    else
        return imageToArray(image);
END_RCPP
}

RcppExport SEXP writeNifti (SEXP _image, SEXP _file)
{
BEGIN_RCPP
    NiftiImage image = retrieveImage(_image);
    
    const int status = nifti_set_filenames(image, as<std::string>(_file).c_str(), false, true);
    if (status != 0)
        throw std::runtime_error("Failed to set filenames for NIfTI object");
    
    nifti_image_write(image);
    
    return R_NilValue;
END_RCPP
}

RcppExport SEXP getXform (SEXP _image, SEXP _preferQuaternion)
{
BEGIN_RCPP
    const NiftiImage image = retrieveImage(_image);
    const bool preferQuaternion = as<bool>(_preferQuaternion);
    
    AffineMatrix matrix;
    
    // No qform or sform so return RAS matrix (NB: other software may assume differently)
    if (image->qform_code <= 0 && image->sform_code <= 0)
        matrix(0,0) = matrix(1,1) = matrix(2,2) = matrix(3,3) = 1.0;
    else if ((preferQuaternion && image->qform_code > 0) || image->sform_code <= 0)
        matrix = AffineMatrix(image->qto_xyz, false);
    else
        matrix = AffineMatrix(image->sto_xyz, false);
    
    return matrix;
END_RCPP
}

RcppExport SEXP regLinear (SEXP _source, SEXP _target, SEXP _type, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _useBlockPercentage, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _init, SEXP _verbose, SEXP _estimateOnly, SEXP _sequentialInit)
{
BEGIN_RCPP
    NiftiImage sourceImage = retrieveImage(_source);
    NiftiImage targetImage = retrieveImage(_target);
    NiftiImage sourceMask = retrieveImage(_sourceMask);
    NiftiImage targetMask = retrieveImage(_targetMask);
    
    if (sourceImage.isNull())
        throw std::runtime_error("Cannot read or retrieve source image");
    if (targetImage.isNull())
        throw std::runtime_error("Cannot read or retrieve target image");
    
    const LinearTransformScope scope = (as<int>(_type) == TYPE_AFFINE ? AffineScope : RigidScope);
    const int interpolation = as<int>(_interpolation);
    const bool symmetric = as<bool>(_symmetric);
    const bool estimateOnly = as<bool>(_estimateOnly);
    const bool sequentialInit = as<bool>(_sequentialInit);
    
    List init(_init);
    List returnValue = List::create(Named("source")=imageToPointer(sourceImage,"Source image"), Named("target")=imageToPointer(targetImage,"Target image"));
    
    if (sourceImage.nDims() == targetImage.nDims())
    {
        AffineMatrix initAffine;
        if (!Rf_isNull(init[0]))
            initAffine = AffineMatrix(SEXP(init[0]));
        else
            initAffine = AffineMatrix(sourceImage, targetImage);
    
        AladinResult result = regAladin(sourceImage, targetImage, scope, symmetric, as<int>(_nLevels), as<int>(_maxIterations), as<int>(_useBlockPercentage), as<int>(_interpolation), sourceMask, targetMask, initAffine, as<bool>(_verbose), estimateOnly);
        
        returnValue["image"] = imageToArray(result.image);
        returnValue["forwardTransforms"] = List::create(result.forwardTransform);
        if (symmetric)
            returnValue["reverseTransforms"] = List::create(result.reverseTransform);
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = List::create(result.iterations);
        
        return returnValue;
    }
    else if (sourceImage.nDims() - targetImage.nDims() == 1)
    {
        const int nReps = sourceImage->dim[sourceImage.nDims()];
        List forwardTransforms(nReps), reverseTransforms(nReps), iterations(nReps);
        NiftiImage finalImage = allocateMultiregResult(sourceImage, targetImage, interpolation != 0);
        AladinResult result;
        for (int i=0; i<nReps; i++)
        {
            NiftiImage currentSource;
            if (sourceImage.nDims() == 3)
                currentSource = sourceImage.slice(i);
            else
                currentSource = sourceImage.volume(i);
            
            AffineMatrix initAffine;
            if (!Rf_isNull(init[i]))
                initAffine = AffineMatrix(SEXP(init[i]));
            else if (sequentialInit && i>0 && result.forwardTransform.isValid())
                initAffine = result.forwardTransform;
            else
                initAffine = AffineMatrix(currentSource, targetImage);
            
            result = regAladin(currentSource, targetImage, scope, symmetric, as<int>(_nLevels), as<int>(_maxIterations), as<int>(_useBlockPercentage), interpolation, sourceMask, targetMask, initAffine, as<bool>(_verbose), estimateOnly);
            
            if (sourceImage.nDims() == 3)
                finalImage.slice(i) = result.image;
            else
                finalImage.volume(i) = result.image;
            
            forwardTransforms[i] = result.forwardTransform;
            if (symmetric)
                reverseTransforms[i] = result.reverseTransform;
            iterations[i] = result.iterations;
        }
        
        returnValue["image"] = imageToArray(finalImage);
        returnValue["forwardTransforms"] = forwardTransforms;
        if (symmetric)
            returnValue["reverseTransforms"] = reverseTransforms;
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = iterations;
        
        return returnValue;
    }
    else
    {
        std::ostringstream message;
        message << "Cannot register a " << sourceImage.nDims() << "D source image to a " << targetImage.nDims() << "D target";
        throw std::runtime_error(message.str());
    }
    
    return R_NilValue;
END_RCPP
}

RcppExport SEXP regNonlinear (SEXP _source, SEXP _target, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _init, SEXP _nBins, SEXP _spacing, SEXP _bendingEnergyWeight, SEXP _jacobianWeight, SEXP _inverseConsistencyWeight, SEXP _verbose, SEXP _estimateOnly, SEXP _sequentialInit)
{
BEGIN_RCPP
    NiftiImage sourceImage = retrieveImage(_source);
    NiftiImage targetImage = retrieveImage(_target);
    NiftiImage sourceMask = retrieveImage(_sourceMask);
    NiftiImage targetMask = retrieveImage(_targetMask);
    
    if (sourceImage.isNull())
        throw std::runtime_error("Cannot read or retrieve source image");
    if (targetImage.isNull())
        throw std::runtime_error("Cannot read or retrieve target image");
    
    const int interpolation = as<int>(_interpolation);
    const bool symmetric = as<bool>(_symmetric);
    const bool estimateOnly = as<bool>(_estimateOnly);
    const bool sequentialInit = as<bool>(_sequentialInit);
    
    List init(_init);
    List returnValue = List::create(Named("source")=imageToPointer(sourceImage,"Source image"), Named("target")=imageToPointer(targetImage,"Target image"));
    
    if (sourceImage.nDims() == targetImage.nDims())
    {
        AffineMatrix initAffine;
        NiftiImage initControl;
        if (!Rf_isNull(init[0]))
        {
            // NB: R code must set the class of an affine appropriately
            RObject initObject(init[0]);
            if (initObject.hasAttribute("class") && as<std::string>(initObject.attr("class")) == "affine")
                initAffine = AffineMatrix(SEXP(initObject));
            else
                initControl = retrieveImage(init[0]);
        }
        else
            initAffine = AffineMatrix(sourceImage, targetImage);
    
        F3dResult result = regF3d(sourceImage, targetImage, as<int>(_nLevels), as<int>(_maxIterations), interpolation, sourceMask, targetMask, initControl, initAffine, as<int>(_nBins), as<float_vector>(_spacing), as<float>(_bendingEnergyWeight), as<float>(_jacobianWeight), as<float>(_inverseConsistencyWeight), symmetric, as<bool>(_verbose), estimateOnly);
        
        returnValue["image"] = imageToArray(result.image);
        returnValue["forwardTransform"] = List::create(imageToPointer(result.forwardTransform, "F3D control points"));
        if (symmetric)
            returnValue["reverseTransforms"] = List::create(imageToPointer(result.reverseTransform, "F3D control points"));
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = List::create(result.iterations);
        
        return returnValue;
    }
    else if (sourceImage.nDims() - targetImage.nDims() == 1)
    {
        const int nReps = sourceImage->dim[sourceImage.nDims()];
        List forwardTransforms(nReps), reverseTransforms(nReps), iterations(nReps);
        NiftiImage finalImage = allocateMultiregResult(sourceImage, targetImage, interpolation != 0);
        F3dResult result;
        for (int i=0; i<nReps; i++)
        {
            NiftiImage currentSource;
            if (sourceImage.nDims() == 3)
                currentSource = sourceImage.slice(i);
            else
                currentSource = sourceImage.volume(i);
            
            AffineMatrix initAffine;
            NiftiImage initControl;
            if (!Rf_isNull(init[i]))
            {
                // NB: R code must set the class of an affine appropriately
                RObject initObject(init[i]);
                if (initObject.hasAttribute("class") && as<std::string>(initObject.attr("class")) == "affine")
                    initAffine = AffineMatrix(SEXP(initObject));
                else
                    initControl = retrieveImage(init[i]);
            }
            else if (sequentialInit && i>0 && !result.forwardTransform.isNull())
                initControl = result.forwardTransform;
            else
                initAffine = AffineMatrix(currentSource, targetImage);
            
            result = regF3d(currentSource, targetImage, as<int>(_nLevels), as<int>(_maxIterations), interpolation, sourceMask, targetMask, initControl, initAffine, as<int>(_nBins), as<float_vector>(_spacing), as<float>(_bendingEnergyWeight), as<float>(_jacobianWeight), as<float>(_inverseConsistencyWeight), symmetric, as<bool>(_verbose), estimateOnly);
            
            if (sourceImage.nDims() == 3)
                finalImage.slice(i) = result.image;
            else
                finalImage.volume(i) = result.image;
            
            forwardTransforms[i] = imageToPointer(result.forwardTransform, "F3D control points");
            if (symmetric)
                reverseTransforms[i] = imageToPointer(result.reverseTransform, "F3D control points");
            iterations[i] = result.iterations;
        }
        
        returnValue["image"] = imageToArray(finalImage);
        returnValue["forwardTransforms"] = forwardTransforms;
        if (symmetric)
            returnValue["reverseTransforms"] = reverseTransforms;
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = iterations;
        
        return returnValue;
    }
    else
    {
        std::ostringstream message;
        message << "Cannot register a " << sourceImage.nDims() << "D source image to a " << targetImage.nDims() << "D target";
        throw std::runtime_error(message.str());
    }
    
    return R_NilValue;
END_RCPP
}

RcppExport SEXP getDeformationField (SEXP _transform)
{
BEGIN_RCPP
    RObject transform(_transform);
    NiftiImage targetImage = retrieveImage(transform.attr("target"));
    
    if (transform.hasAttribute("class") && as<std::string>(transform.attr("class")) == "affine")
    {
        AffineMatrix affine = AffineMatrix(SEXP(transform));
        DeformationField field(targetImage, affine);
        return (imageToPointer(field.getFieldImage(), "Deformation field"));
    }
    else
    {
        NiftiImage transformationImage = retrieveImage(_transform);
        DeformationField field(targetImage, transformationImage);
        return (imageToPointer(field.getFieldImage(), "Deformation field"));
    }
    
    return R_NilValue;
END_RCPP
}

RcppExport SEXP transformPoints (SEXP _transform, SEXP _points, SEXP _nearest)
{
BEGIN_RCPP
    NiftiImage transformationImage = retrieveImage(_transform);
    RObject transform(_transform);
    NiftiImage targetImage = retrieveImage(transform.attr("target"));
    DeformationField deformationField(targetImage, transformationImage);
    NumericMatrix points(_points);
    List result(points.nrow());
    const bool nearest = as<bool>(_nearest);
    
    if (points.ncol() == 2)
    {
        for (int i=0; i<points.nrow(); i++)
        {
            Eigen::Vector2d point;
            point[0] = points(i, 0);
            point[1] = points(i, 1);
            result[i] = deformationField.findPoint(point, nearest);
        }
    }
    else if (points.ncol() == 3)
    {
        for (int i=0; i<points.nrow(); i++)
        {
            Eigen::Vector3d point;
            point[0] = points(i, 0);
            point[1] = points(i, 1);
            point[2] = points(i, 2);
            result[i] = deformationField.findPoint(point, nearest);
        }
    }
    else
        throw std::runtime_error("Points matrix should have 2 or 3 columns");
    
    return result;
END_RCPP
}
