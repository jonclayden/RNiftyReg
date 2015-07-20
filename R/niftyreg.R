niftyreg <- function (source, target, scope = c("affine","rigid","nonlinear"), init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, estimateOnly = FALSE, ...)
{
    if (missing(source) || missing(target))
        stop("Source and target images must be given")
    
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
        stop("Source and target images must be given")
    if (!(nSourceDim %in% 2:4))
        stop("Only 2D, 3D or 4D source images may be used")
    if (!(nTargetDim %in% 2:3))
        stop("Only 2D or 3D target images may be used")
    if (nSourceDim == 4 && nTargetDim == 2)
        stop("4D to 2D registration cannot be performed")
    if (any(dim(source) < 4) || any(dim(target) < 4))
        stop("Images of fewer than 4 voxels in any dimension cannot be registered")
}

niftyreg.linear <- function (source, target, scope = c("affine","rigid"), init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, nLevels = 3L, maxIterations = 5L, useBlockPercentage = 50L, interpolation = 3L, verbose = FALSE, estimateOnly = FALSE, sequentialInit = FALSE)
{
    nSourceDim <- length(dim(source))
    nTargetDim <- length(dim(target))
    .checkImages(source, target)
    
    if (!(interpolation %in% c(0,1,3)))
        stop("Final interpolation specifier must be 0, 1 or 3")
    
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
            stop("Linear registration can only be initialised with an affine matrix")
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
        stop("Penalty term weights must be nonnegative")
    if (bendingEnergyWeight + jacobianWeight > 1)
        stop("Penalty term weights cannot add up to more than 1")
    if (!(interpolation %in% c(0,1,3)))
        stop("Final interpolation specifier must be 0, 1 or 3")
    
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
                stop("Initial control point images must all use the same grid")
            finalSpacing <<- currentSpacing
            spacingUnit <<- "mm"
            return (x)
        }
        else if (!is.matrix(x) || !isTRUE(all.equal(dim(x), c(4,4))))
            stop("Initial transform should be a control point image or affine matrix")
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
        attr(result, "source") <- object$source[[i]]
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
        result <- object$reverseTransforms[[i]]
        attr(result, "source") <- object$target
        attr(result, "target") <- object$source[[i]]
        return (result)
    }
}
