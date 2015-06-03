#include <RcppEigen.h>

#include "config.h"
#include "nifti_image.h"
#include "aladin.h"

// Registration types (degrees of freedom)
#define TYPE_RIGID  0
#define TYPE_AFFINE 1

using namespace Rcpp;

RcppExport SEXP regLinear (SEXP _source, SEXP _target, SEXP _type, SEXP _symmetric, SEXP _nLevels, SEXP _maxIterations, SEXP _useBlockPercentage, SEXP _interpolation, SEXP _sourceMask, SEXP _targetMask, SEXP _initAffine, SEXP _verbose, SEXP _estimateOnly)
{
BEGIN_RCPP
    nifti_image *sourceImage = retrieveImage(_source);
    nifti_image *targetImage = retrieveImage(_target);
    nifti_image *sourceMask = retrieveImage(_sourceMask);
    nifti_image *targetMask = retrieveImage(_targetMask);
    
    LinearTransformScope scope = (as<int>(_type) == TYPE_AFFINE ? AffineScope : RigidScope);
    
    AffineMatrix *initAffine = NULL;
    if (!Rf_isNull(_initAffine))
        initAffine = new AffineMatrix(_initAffine);
    
    AladinResult result = regAladin(sourceImage, targetImage, scope, as<bool>(_symmetric), as<int>(_nLevels), as<int>(_maxIterations), as<int>(_useBlockPercentage), as<int>(_interpolation), sourceMask, targetMask, initAffine, as<bool>(_verbose), as<bool>(_estimateOnly));
    
    nifti_image_free(sourceImage);
    nifti_image_free(targetImage);
    nifti_image_free(sourceMask);
    nifti_image_free(targetMask);
    
    delete initAffine;
    
    return R_NilValue;
END_RCPP
}
