.fixTypes <- function (image)
{
    if (is.null(image))
        return (NULL)
    else
    {
        integerSlots <- c("dim_", "intent_code", "datatype", "bitpix", "slice_start", "slice_end", "slice_code", "xyzt_units", "qform_code", "sform_code")
        numericSlots <- c("intent_p1", "intent_p2", "intent_p3", "slice_duration", "pixdim", "vox_offset", "scl_slope", "scl_inter", "toffset", "cal_max", "cal_min", "quatern_b", "quatern_c", "quatern_d", "qoffset_x", "qoffset_y", "qoffset_z", "srow_x", "srow_y", "srow_z")
        
        for (i in seq_along(integerSlots))
            slot(image, integerSlots[i]) <- as.integer(slot(image, integerSlots[i]))
        
        for (i in seq_along(numericSlots))
            slot(image, numericSlots[i]) <- as.numeric(slot(image, numericSlots[i]))
        
        invalidDatatypes <- c("COMPLEX64", "RGB24", "INT64", "UINT64", "FLOAT128", "COMPLEX128", "COMPLEX256", "RGBA32")
        doubleDatatypes <- c("FLOAT32", "FLOAT64")
        
        datatypeName <- convert.datatype(image@datatype)
        if (datatypeName %in% invalidDatatypes)
            report(OL$Error, "RNiftyReg does not support the \"", datatypeName, "\" image data type")
        else if (datatypeName %in% doubleDatatypes)
            storage.mode(image@.Data) <- "double"
        else
            storage.mode(image@.Data) <- "integer"
        
        return (image)
    }
}

.createControlPointImage <- function (data, dim, pixdim, xform)
{
    dim(data) <- dim
    extendedDim <- c(length(dim), dim, rep(1,7-length(dim)))
    
    # The validity method for the "nifti" class requires this
    if (xform[1] != 0 && xform[9] == 0)
        xform[9] <- 1
    
    image <- new("nifti", .Data=data, dim_=extendedDim, datatype=64L, bitpix=64, pixdim=c(xform[9],pixdim,1,1,0,0), xyzt_units=0, qform_code=xform[1], sform_code=xform[2], quatern_b=xform[3], quatern_c=xform[4], quatern_d=xform[5], qoffset_x=xform[6], qoffset_y=xform[7], qoffset_z=xform[8], srow_x=xform[10:13], srow_y=xform[14:17], srow_z=xform[18:21], cal_min=min(data,na.rm=TRUE), cal_max=max(data,na.rm=TRUE))
    
    return (image)
}

.setImageMetadata <- function (image, referenceImage, finalInterpolation)
{
    if (is.null(image))
        return (NULL)
    
    image@cal_min <- min(image@.Data, na.rm=TRUE)
    image@cal_max <- max(image@.Data, na.rm=TRUE)
    image@scl_slope <- referenceImage@scl_slope
    image@scl_inter <- referenceImage@scl_inter
    
    if (finalInterpolation == 0)
    {
        image@datatype <- referenceImage@datatype
        image@bitpix <- as.numeric(referenceImage@bitpix)
    }
    else
    {
        image@datatype <- 64L
        image@bitpix <- 64
    }
    image@data_type <- convert.datatype(image@datatype)
    
    return (image)
}

pixdim <- function (object)
{
    if (!is.null(attr(object, "pixdim")))
        return (attr(object, "pixdim"))
    else if (!is.null(dim(object)))
        return (rep(1, length(dim(object))))
    else
        return (1)
}

niftyreg <- function (source, target, scope = c("affine","rigid","nonlinear"), init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, estimateOnly = FALSE, ...)
{
    if (missing(source) || missing(target))
        report(OL$Error, "Source and target images must be given")
    
    scope <- match.arg(scope)
    if (scope == "nonlinear")
        niftyreg.nonlinear(source, target, init, sourceMask, targetMask, symmetric=symmetric, estimateOnly=estimateOnly, ...)
    else
        niftyreg.linear(source, target, scope, init, sourceMask, targetMask, symmetric=symmetric, estimateOnly=estimateOnly, ...)
}

# Standard set of preregistration checks for source and target images
.checkImages <- function (source, target)
{
    nSourceDim <- length(dim(source))
    nTargetDim <- length(dim(target))
    
    if (missing(source) || missing(target))
        report(OL$Error, "Source and target images must be given")
    if (!(nSourceDim %in% 2:4))
        report(OL$Error, "Only 2D, 3D or 4D source images may be used")
    if (!(nTargetDim %in% 2:3))
        report(OL$Error, "Only 2D or 3D target images may be used")
    if (nSourceDim == 4 && nTargetDim == 2)
        report(OL$Error, "4D to 2D registration cannot be performed")
    if (any(dim(source) < 4) || any(dim(target) < 4))
        report(OL$Error, "Images of fewer than 4 voxels in any dimension cannot be registered")
}

niftyreg.linear <- function (source, target, scope = c("affine","rigid"), init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, nLevels = 3L, maxIterations = 5L, useBlockPercentage = 50L, interpolation = 3L, verbose = FALSE, estimateOnly = FALSE, sequentialInit = FALSE)
{
    nSourceDim <- length(dim(source))
    nTargetDim <- length(dim(target))
    .checkImages(source, target)
    
    if (!(interpolation %in% c(0,1,3)))
        report(OL$Error, "Final interpolation specifier must be 0, 1 or 3")
    
    scope <- match.arg(scope)
    nReps <- ifelse(nSourceDim > nTargetDim, dim(source)[nSourceDim], 1L)
    
    if (!is.list(init))
        init <- list(init)
    if (length(init) != nReps)
    {
        if (sequentialInit)
            init <- c(init, rep(list(NULL),nReps-length(init)))
        else
            init <- rep(init, length.out=nReps)
    }
    init <- lapply(init, function(x) {
        if (is.null(x))
            return (x)
        else if (!is.matrix(x) || !isTRUE(all.equal(dim(x), c(4,4))))
            report(OL$Error, "Linear registration can only be initialised with an affine matrix")
        else if (!is.null(attr(x,"affineType")) && attr(x,"affineType") != "niftyreg")
            return (convertAffine(x, source, target, "niftyreg"))
    })
    
    result <- .Call("regLinear", source, target, ifelse(scope=="affine",1L,0L), symmetric, nLevels, maxIterations, useBlockPercentage, interpolation, sourceMask, targetMask, init, verbose, estimateOnly, sequentialInit, PACKAGE="RNiftyReg")
    class(result) <- "niftyreg"
    
    return (result)
}

niftyreg.nonlinear <- function (source, target, init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, nLevels = 3L, maxIterations = 300L, nBins = 64L, bendingEnergyWeight = 0.005, jacobianWeight = 0, inverseConsistencyWeight = 0.01, finalSpacing = c(5,5,5), spacingUnit = c("vox","mm"), interpolation = 3L, verbose = FALSE, estimateOnly = FALSE, sequentialInit = FALSE)
{
    nSourceDim <- length(dim(source))
    nTargetDim <- length(dim(target))
    .checkImages(source, target)
    
    if (any(c(bendingEnergyWeight,jacobianWeight,inverseConsistencyWeight) < 0))
        report(OL$Error, "Penalty term weights must be nonnegative")
    if (bendingEnergyWeight + jacobianWeight > 1)
        report(OL$Error, "Penalty term weights cannot add up to more than 1")
    if (!(interpolation %in% c(0,1,3)))
        report(OL$Error, "Final interpolation specifier must be 0, 1 or 3")
    
    if (nLevels == 0)
        symmetric <- FALSE
    
    nReps <- ifelse(nSourceDim > nTargetDim, dim(source)[nSourceDim], 1L)
    spacingUnit <- match.arg(spacingUnit)
    spacingChanged <- FALSE
    
    if (!is.list(init))
        init <- list(init)
    if (length(init) != nReps)
    {
        if (sequentialInit)
            init <- c(init, rep(list(NULL),nReps-length(init)))
        else
            init <- rep(init, length.out=nReps)
    }
    init <- lapply(init, function(x) {
        if (is.null(x))
            return (x)
        else if (!is.null(attr(x, ".nifti_image_ptr")))
        {
            currentSpacing <- pixdim(x)[1:3] / 2^max(0,nLevels-1)
            if (spacingChanged && !isTRUE(all.equal(currentSpacing, finalSpacing)))
                report(OL$Error, "Initial control point images must all use the same grid")
            finalSpacing <<- currentSpacing
            spacingUnit <<- "mm"
            return (x)
        }
        else if (!is.matrix(x) || !isTRUE(all.equal(dim(x), c(4,4))))
            report(OL$Error, "Initial transform should be a control point image or affine matrix")
        else if (!is.null(attr(x,"affineType")) && attr(x,"affineType") != "niftyreg")
            return (convertAffine(x, source, target, "niftyreg"))
    })
    
    if (spacingUnit == "vox")
        finalSpacing <- finalSpacing * abs(pixdim(target)[1:min(3,nTargetDim)])
    
    if (nTargetDim == 2)
        finalSpacing <- c(finalSpacing[1:2], 1)
    else
        finalSpacing <- finalSpacing[1:3]
    
    result <- .Call("regNonlinear", source, target, symmetric, nLevels, maxIterations, interpolation, sourceMask, targetMask, init, nBins, finalSpacing, bendingEnergyWeight, jacobianWeight, inverseConsistencyWeight, verbose, estimateOnly, sequentialInit, PACKAGE="RNiftyReg")
    class(result) <- "niftyreg"
    
    return (result)
}

forward <- function (object, ...)
{
    UseMethod("forward")
}

forward.niftyreg <- function (object, i = 1, ...)
{
    if (is.null(object$forwardTransforms))
        return (NULL)
    else
    {
        result <- object$forwardTransforms[[i]]
        attr(result, "source") <- object$source
        attr(result, "target") <- object$target
        return (result)
    }
}

reverse <- function (object, ...)
{
    UseMethod("reverse")
}

reverse.niftyreg <- function (object, i = 1, ...)
{
    if (is.null(object$reverseTransforms))
        return (NULL)
    else
    {
        result <- object$forwardTransforms[[i]]
        attr(result, "source") <- object$source
        attr(result, "target") <- object$target
        return (result)
    }
}

applyAffine <- function (affine, source, target, affineType = NULL, finalInterpolation = 3)
{
    if (!is.matrix(affine) || !isTRUE(all.equal(dim(affine), c(4,4))))
        report(OL$Error, "Specified affine matrix is not valid")
    
    if (is.null(affineType))
    {
        affineType <- attr(affine, "affineType")
        if (is.null(affineType))
            report(OL$Error, "The current affine type was not specified and is not stored with the matrix")
    }
    else
        attr(affine, "affineType") <- affineType
    
    return (niftyreg.linear(source, target, targetMask=NULL, initAffine=affine, scope="affine", nLevels=0, finalInterpolation=finalInterpolation, verbose=FALSE, estimateOnly=FALSE))
}

applyControlPoints <- function (controlPointImage, source, target, finalInterpolation = 3)
{
    return (niftyreg.nonlinear(source, target, targetMask=NULL, initControl=controlPointImage, symmetric=FALSE, nLevels=0, finalInterpolation=finalInterpolation, verbose=FALSE, estimateOnly=FALSE))
}

getDeformationField <- function (target, affine = NULL, controlPointImage = NULL, jacobian = TRUE)
{
    if (missing(target))
        report(OL$Error, "Target image must be given")
    else
        target <- as(target, "nifti")
    if (is.null(affine) && is.null(controlPointImage))
        report(OL$Error, "Affine matrix or control point image must be specified")
    if (!is.null(controlPointImage))
        controlPointImage <- as(controlPointImage, "nifti")
    
    returnValue <- .Call("get_deformation_field_R", affine, .fixTypes(controlPointImage), .fixTypes(target), jacobian, PACKAGE="RNiftyReg")
    
    nDims <- target@dim_[1]
    dimIndex <- 1 + seq_len(nDims)
    padding <- rep(1, 4-nDims)
    
    result <- list(deformationField=.createControlPointImage(returnValue[[1]], c(target@dim_[dimIndex],padding,nDims), c(target@pixdim[dimIndex],padding,1), returnValue[[2]]))
    if (jacobian)
        result$jacobian <- .createControlPointImage(returnValue[[3]], target@dim_[dimIndex], target@pixdim[dimIndex], returnValue[[4]])
    
    return (result)
}
