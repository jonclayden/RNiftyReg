#include <RcppEigen.h>

#include "config.h"
#include "NiftiImage.h"
#include "aladin.h"

// Registration types (degrees of freedom)
#define TYPE_RIGID  0
#define TYPE_AFFINE 1

using namespace Rcpp;

RcppExport SEXP regLinear (SEXP _source, SEXP _target, SEXP _type, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _useBlockPercentage, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _initAffine, SEXP _verbose, SEXP _estimateOnly)
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
    
    LinearTransformScope scope = (as<int>(_type) == TYPE_AFFINE ? AffineScope : RigidScope);
    
    AffineMatrix initAffine;
    if (Rf_isNull(_initAffine))
        initAffine = AffineMatrix(sourceImage, targetImage);
    else
        initAffine = AffineMatrix(_initAffine);
    
    AladinResult result = regAladin(sourceImage, targetImage, scope, as<bool>(_symmetric), as<int>(_nLevels), as<int>(_maxIterations), as<int>(_useBlockPercentage), as<int>(_interpolation), sourceMask, targetMask, initAffine, as<bool>(_verbose), as<bool>(_estimateOnly));
    
    return List::create(Named("image")=imageToArray(result.image), Named("affine")=result.affine, Named("iterations")=result.iterations);
END_RCPP
}
