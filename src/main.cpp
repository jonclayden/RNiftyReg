#include <RcppEigen.h>

#include "config.h"
#include "NiftiImage.h"
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
    
    if (sourceImage.nDims() == targetImage.nDims())
    {
        AffineMatrix initAffine;
        if (!Rf_isNull(init[0]))
            initAffine = AffineMatrix(SEXP(init[0]));
        else
            initAffine = AffineMatrix(sourceImage, targetImage);
    
        AladinResult result = regAladin(sourceImage, targetImage, scope, symmetric, as<int>(_nLevels), as<int>(_maxIterations), as<int>(_useBlockPercentage), as<int>(_interpolation), sourceMask, targetMask, initAffine, as<bool>(_verbose), estimateOnly);
        
        List returnValue = List::create(Named("image")=imageToArray(result.image), Named("forwardTransforms")=List::create(result.forwardTransform), Named("iterations")=List::create(result.iterations));
        
        if (symmetric)
            returnValue["reverseTransforms"] = List::create(result.reverseTransform);
        else
            returnValue["reverseTransforms"] = R_NilValue;
        
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
            reverseTransforms[i] = result.reverseTransform;
            iterations[i] = result.iterations;
        }
        
        List returnValue = List::create(Named("image")=imageToArray(finalImage), Named("forwardTransforms")=forwardTransforms, Named("iterations")=iterations);
        
        if (symmetric)
            returnValue["reverseTransforms"] = reverseTransforms;
        else
            returnValue["reverseTransforms"] = R_NilValue;
        
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

RcppExport SEXP regNonlinear (SEXP _source, SEXP _target, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _useBlockPercentage, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _initControlPoints, SEXP _initAffine, SEXP _nBins, SEXP _spacing, SEXP _bendingEnergyWeight, SEXP _jacobianWeight, SEXP _inverseConsistencyWeight, SEXP _verbose, SEXP _estimateOnly)
{
BEGIN_RCPP
    NiftiImage sourceImage = retrieveImage(_source);
    NiftiImage targetImage = retrieveImage(_target);
    NiftiImage sourceMask = retrieveImage(_sourceMask);
    NiftiImage targetMask = retrieveImage(_targetMask);
    NiftiImage initControlPoints = retrieveImage(_initControlPoints);
    
    if (sourceImage.isNull())
        throw std::runtime_error("Cannot read or retrieve source image");
    if (targetImage.isNull())
        throw std::runtime_error("Cannot read or retrieve target image");
    
    AffineMatrix initAffine;
    if (!Rf_isNull(_initAffine))
        initAffine = AffineMatrix(_initAffine);
    else if (initControlPoints.isNull())
        initAffine = AffineMatrix(sourceImage, targetImage);
    
    F3dResult result = regF3d(sourceImage, targetImage, as<int>(_nLevels), as<int>(_maxIterations), as<int>(_interpolation), sourceMask, targetMask, initControlPoints, initAffine, as<int>(_nBins), as<float_vector>(_spacing), as<float>(_bendingEnergyWeight), as<float>(_jacobianWeight), as<float>(_inverseConsistencyWeight), as<bool>(_symmetric), as<bool>(_verbose), as<bool>(_estimateOnly));
    
    return List::create(Named("image")=imageToArray(result.forwardImage), Named("control")=imageToPointer(result.forwardControlPoints, "f3d control points"), Named("reverseImage")=imageToArray(result.reverseImage), Named("reverseControl")=imageToPointer(result.reverseControlPoints, "f3d control points"), Named("iterations")=result.iterations);
END_RCPP
}
