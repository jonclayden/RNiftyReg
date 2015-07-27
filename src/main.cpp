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

RcppExport SEXP retrieveImage (SEXP _image)
{
BEGIN_RCPP
    NiftiImage image(_image);
    return image.toPointer("NIfTI image");
END_RCPP
}

RcppExport SEXP readNifti (SEXP _file, SEXP _internal)
{
BEGIN_RCPP
    NiftiImage image(_file);
    
    if (as<bool>(_internal))
        return image.toPointer("NIfTI image");
    else
        return image.toArray();
END_RCPP
}

RcppExport SEXP writeNifti (SEXP _image, SEXP _file)
{
BEGIN_RCPP
    NiftiImage image(_image);
    
    const int status = nifti_set_filenames(image, as<std::string>(_file).c_str(), false, true);
    if (status != 0)
        throw std::runtime_error("Failed to set filenames for NIfTI object");
    
    nifti_image_write(image);
    
    return R_NilValue;
END_RCPP
}

RcppExport SEXP updateNifti (SEXP _image, SEXP _reference)
{
BEGIN_RCPP
    const NiftiImage reference(_reference);
    RObject object(_image);
    
    if (!reference.isNull())
    {
        NiftiImage *updatedImage = new NiftiImage(reference);
        updatedImage->update(_image);
        // NiftiImage *updatedImage = new NiftiImage(reference, _image);
        updatedImage->setPersistence(true);
        XPtr<NiftiImage> xptr(updatedImage);
        R_RegisterCFinalizerEx(SEXP(xptr), &finaliseNiftiImage, FALSE);
        object.attr(".nifti_image_ptr") = xptr;
    }
    
    return object;
END_RCPP
}

RcppExport SEXP dumpNifti (SEXP _image)
{
BEGIN_RCPP
    const NiftiImage image(_image);
    return image.headerToList();
END_RCPP
}

RcppExport SEXP getXform (SEXP _image, SEXP _preferQuaternion)
{
BEGIN_RCPP
    const NiftiImage image(_image);
    const bool preferQuaternion = as<bool>(_preferQuaternion);
    
    AffineMatrix matrix(image.xform(preferQuaternion), false);
    return matrix;
END_RCPP
}

void checkImages (const NiftiImage &sourceImage, const NiftiImage &targetImage)
{
    if (sourceImage.isNull())
        throw std::runtime_error("Cannot read or retrieve source image");
    if (targetImage.isNull())
        throw std::runtime_error("Cannot read or retrieve target image");
    
    const int nSourceDim = sourceImage.nDims();
    const int nTargetDim = targetImage.nDims();
    
    if (nSourceDim < 2 || nSourceDim > 4)
        throw std::runtime_error("Source image should have 2, 3 or 4 dimensions");
    if (nTargetDim < 2 || nTargetDim > 3)
        throw std::runtime_error("Target image should have 2 or 3 dimensions");
}

RcppExport SEXP regLinear (SEXP _source, SEXP _target, SEXP _type, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _useBlockPercentage, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _init, SEXP _verbose, SEXP _estimateOnly, SEXP _sequentialInit)
{
BEGIN_RCPP
    NiftiImage sourceImage(_source);
    NiftiImage targetImage(_target);
    NiftiImage sourceMask(_sourceMask);
    NiftiImage targetMask(_targetMask);
    
    checkImages(sourceImage.drop(), targetImage.drop());
    
    const LinearTransformScope scope = (as<int>(_type) == TYPE_AFFINE ? AffineScope : RigidScope);
    const int interpolation = as<int>(_interpolation);
    const bool symmetric = as<bool>(_symmetric);
    const bool estimateOnly = as<bool>(_estimateOnly);
    const bool sequentialInit = as<bool>(_sequentialInit);
    
    List init(_init);
    List returnValue;
    
    if (sourceImage.nDims() == targetImage.nDims())
    {
        AffineMatrix initAffine;
        if (!Rf_isNull(init[0]))
            initAffine = AffineMatrix(SEXP(init[0]));
        else
            initAffine = AffineMatrix(sourceImage, targetImage);
    
        AladinResult result = regAladin(sourceImage, targetImage, scope, symmetric, as<int>(_nLevels), as<int>(_maxIterations), as<int>(_useBlockPercentage), as<int>(_interpolation), sourceMask, targetMask, initAffine, as<bool>(_verbose), estimateOnly);
        
        returnValue["image"] = result.image.toArray();
        returnValue["forwardTransforms"] = List::create(result.forwardTransform);
        if (symmetric)
            returnValue["reverseTransforms"] = List::create(result.reverseTransform);
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = List::create(result.iterations);
        returnValue["source"] = List::create(sourceImage.toPointer("Source image"));
        returnValue["target"] = targetImage.toPointer("Target image");
        
        return returnValue;
    }
    else if (sourceImage.nDims() - targetImage.nDims() == 1)
    {
        const int nReps = sourceImage->dim[sourceImage.nDims()];
        List forwardTransforms(nReps), reverseTransforms(nReps), iterations(nReps), sourceImages(nReps);
        NiftiImage finalImage = allocateMultiregResult(sourceImage, targetImage, interpolation != 0);
        AladinResult result;
        for (int i=0; i<nReps; i++)
        {
            NiftiImage currentSource;
            if (sourceImage.nDims() == 3)
                currentSource = sourceImage.slice(i);
            else
                currentSource = sourceImage.volume(i);
            sourceImages[i] = currentSource.toPointer("Source image");
            
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
        
        returnValue["image"] = finalImage.toArray();
        returnValue["forwardTransforms"] = forwardTransforms;
        if (symmetric)
            returnValue["reverseTransforms"] = reverseTransforms;
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = iterations;
        returnValue["source"] = sourceImages;
        returnValue["target"] = targetImage.toPointer("Target image");
        
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

RcppExport SEXP regNonlinear (SEXP _source, SEXP _target, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _init, SEXP _nBins, SEXP _spacing, SEXP _bendingEnergyWeight, SEXP _jacobianWeight, SEXP _verbose, SEXP _estimateOnly, SEXP _sequentialInit)
{
BEGIN_RCPP
    NiftiImage sourceImage(_source);
    NiftiImage targetImage(_target);
    NiftiImage sourceMask(_sourceMask);
    NiftiImage targetMask(_targetMask);
    
    checkImages(sourceImage.drop(), targetImage.drop());
    
    const int interpolation = as<int>(_interpolation);
    const bool symmetric = as<bool>(_symmetric);
    const bool estimateOnly = as<bool>(_estimateOnly);
    const bool sequentialInit = as<bool>(_sequentialInit);
    
    List init(_init);
    List returnValue;
    
    if (sourceImage.nDims() == targetImage.nDims())
    {
        AffineMatrix initAffine;
        NiftiImage initControl;
        if (!Rf_isNull(init[0]))
        {
            // NB: R code must set the class of an affine appropriately
            RObject initObject(init[0]);
            if (initObject.inherits("affine"))
                initAffine = AffineMatrix(SEXP(initObject));
            else
                initControl = NiftiImage(SEXP(init[0]));
        }
        else
            initAffine = AffineMatrix(sourceImage, targetImage);
    
        F3dResult result = regF3d(sourceImage, targetImage, as<int>(_nLevels), as<int>(_maxIterations), interpolation, sourceMask, targetMask, initControl, initAffine, as<int>(_nBins), as<float_vector>(_spacing), as<float>(_bendingEnergyWeight), as<float>(_jacobianWeight), symmetric, as<bool>(_verbose), estimateOnly);
        
        returnValue["image"] = result.image.toArray();
        returnValue["forwardTransforms"] = List::create(result.forwardTransform.toPointer("F3D control points"));
        if (symmetric)
            returnValue["reverseTransforms"] = List::create(result.reverseTransform.toPointer("F3D control points"));
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = List::create(result.iterations);
        returnValue["source"] = List::create(sourceImage.toPointer("Source image"));
        returnValue["target"] = targetImage.toPointer("Target image");
        
        return returnValue;
    }
    else if (sourceImage.nDims() - targetImage.nDims() == 1)
    {
        const int nReps = sourceImage->dim[sourceImage.nDims()];
        List forwardTransforms(nReps), reverseTransforms(nReps), iterations(nReps), sourceImages(nReps);
        NiftiImage finalImage = allocateMultiregResult(sourceImage, targetImage, interpolation != 0);
        F3dResult result;
        for (int i=0; i<nReps; i++)
        {
            NiftiImage currentSource;
            if (sourceImage.nDims() == 3)
                currentSource = sourceImage.slice(i);
            else
                currentSource = sourceImage.volume(i);
            sourceImages[i] = currentSource.toPointer("Source image");
            
            AffineMatrix initAffine;
            NiftiImage initControl;
            if (!Rf_isNull(init[i]))
            {
                // NB: R code must set the class of an affine appropriately
                RObject initObject(init[i]);
                if (initObject.inherits("affine"))
                    initAffine = AffineMatrix(SEXP(initObject));
                else
                    initControl = NiftiImage(SEXP(init[i]));
            }
            else if (sequentialInit && i>0 && !result.forwardTransform.isNull())
                initControl = result.forwardTransform;
            else
                initAffine = AffineMatrix(currentSource, targetImage);
            
            result = regF3d(currentSource, targetImage, as<int>(_nLevels), as<int>(_maxIterations), interpolation, sourceMask, targetMask, initControl, initAffine, as<int>(_nBins), as<float_vector>(_spacing), as<float>(_bendingEnergyWeight), as<float>(_jacobianWeight), symmetric, as<bool>(_verbose), estimateOnly);
            
            if (sourceImage.nDims() == 3)
                finalImage.slice(i) = result.image;
            else
                finalImage.volume(i) = result.image;
            
            forwardTransforms[i] = result.forwardTransform.toPointer("F3D control points");
            if (symmetric)
                reverseTransforms[i] = result.reverseTransform.toPointer("F3D control points");
            iterations[i] = result.iterations;
        }
        
        returnValue["image"] = finalImage.toArray();
        returnValue["forwardTransforms"] = forwardTransforms;
        if (symmetric)
            returnValue["reverseTransforms"] = reverseTransforms;
        else
            returnValue["reverseTransforms"] = R_NilValue;
        returnValue["iterations"] = iterations;
        returnValue["source"] = sourceImages;
        returnValue["target"] = targetImage.toPointer("Target image");
        
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

RcppExport SEXP getDeformationField (SEXP _transform, SEXP _jacobian)
{
BEGIN_RCPP
    RObject transform(_transform);
    RObject result;
    NiftiImage targetImage(SEXP(transform.attr("target")));
    DeformationField field;
    
    if (transform.inherits("affine"))
    {
        AffineMatrix affine = AffineMatrix(SEXP(transform));
        field = DeformationField(targetImage, affine);
    }
    else
    {
        NiftiImage transformationImage(_transform);
        field = DeformationField(targetImage, transformationImage);
    }
    
    result = field.getFieldImage().toPointer("Deformation field");
    result.attr("source") = transform.attr("source");
    result.attr("target") = transform.attr("target");
    
    if (as<bool>(_jacobian))
        result.attr("jacobian") = field.getJacobian().toPointer("Jacobian of deformation field");
    
    return result;
END_RCPP
}

RcppExport SEXP transformPoints (SEXP _transform, SEXP _points, SEXP _nearest)
{
BEGIN_RCPP
    NiftiImage transformationImage(_transform);
    RObject transform(_transform);
    NiftiImage sourceImage(SEXP(transform.attr("source")));
    NiftiImage targetImage(SEXP(transform.attr("target")));
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
            result[i] = deformationField.findPoint(sourceImage, point, nearest);
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
            result[i] = deformationField.findPoint(sourceImage, point, nearest);
        }
    }
    else
        throw std::runtime_error("Points matrix should have 2 or 3 columns");
    
    return result;
END_RCPP
}

RcppExport SEXP pointerToArray (SEXP _image)
{
BEGIN_RCPP
    NiftiImage image(_image);
    return image.toArray();
END_RCPP
}

RcppExport SEXP halfTransform (SEXP _transform)
{
BEGIN_RCPP
    RObject transform(_transform);
    RObject result;
    if (transform.inherits("affine"))
    {
        Eigen::MatrixXd matrix = as<Eigen::MatrixXd>(_transform);
        matrix = (matrix.log() * 0.5).exp();
        result = AffineMatrix(matrix);
    }
    else
    {
        NiftiImage transformationImage(_transform);
        switch (reg_round(transformationImage->intent_p1))
        {
            case SPLINE_GRID:
            reg_getDisplacementFromDeformation(transformationImage);
            reg_tools_multiplyValueToImage(transformationImage,transformationImage,0.5f);
            reg_getDeformationFromDisplacement(transformationImage);
            break;
            
            case DEF_FIELD:
            reg_getDisplacementFromDeformation(transformationImage);
            reg_tools_multiplyValueToImage(transformationImage,transformationImage,0.5f);
            reg_getDeformationFromDisplacement(transformationImage);
            break;
            
            case DISP_FIELD:
            reg_tools_multiplyValueToImage(transformationImage,transformationImage,0.5f);
            break;
            
            case SPLINE_VEL_GRID:
            reg_getDisplacementFromDeformation(transformationImage);
            reg_tools_multiplyValueToImage(transformationImage,transformationImage,0.5f);
            reg_getDeformationFromDisplacement(transformationImage);
            --transformationImage->intent_p2;
            if (transformationImage->num_ext>1)
                --transformationImage->num_ext;
            break;
            
            case DEF_VEL_FIELD:
            reg_getDisplacementFromDeformation(transformationImage);
            reg_tools_multiplyValueToImage(transformationImage,transformationImage,0.5f);
            reg_getDeformationFromDisplacement(transformationImage);
            --transformationImage->intent_p2;
            break;
            
            case DISP_VEL_FIELD:
            reg_tools_multiplyValueToImage(transformationImage,transformationImage,0.5f);
            --transformationImage->intent_p2;
            break;
            
            default:
            throw std::runtime_error("The specified transformation image is not valid or not supported");
        }
        
        result = transformationImage.toPointer("F3D transformation");
    }
    
    result.attr("source") = transform.attr("source");
    result.attr("target") = transform.attr("target");
    
    return result;
END_RCPP
}
