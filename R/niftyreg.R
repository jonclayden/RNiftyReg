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
    image <- new("nifti", .Data=data, dim_=extendedDim, datatype=64L, bitpix=64, pixdim=c(xform[9],pixdim,1,1,0,0), xyzt_units=0, qform_code=xform[1], sform_code=xform[2], quatern_b=xform[3], quatern_c=xform[4], quatern_d=xform[5], qoffset_x=xform[6], qoffset_y=xform[7], qoffset_z=xform[8], srow_x=xform[10:13], srow_y=xform[14:17], srow_z=xform[18:21], cal_min=min(data,na.rm=TRUE), cal_max=max(data,na.rm=TRUE))
    
    return (image)
}

.setImageMetadata <- function (image, referenceImage, finalInterpolation)
{
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

niftyreg <- function (source, target, targetMask = NULL, initAffine = NULL, scope = c("affine","rigid","nonlinear"), ...)
{
    if (missing(source) || missing(target))
        report(OL$Error, "Source and target images must be given")
    
    scope <- match.arg(scope)
    if (scope == "nonlinear")
        niftyreg.nonlinear(source, target, targetMask, initAffine, ...)
    else
        niftyreg.linear(source, target, targetMask, initAffine, scope, ...)
}

niftyreg.linear <- function (source, target, targetMask = NULL, initAffine = NULL, scope = c("affine","rigid"), nLevels = 3, maxIterations = 5, useBlockPercentage = 50, finalInterpolation = 3, verbose = FALSE)
{
    if (!require("oro.nifti"))
        report(OL$Error, "The \"oro.nifti\" package is required")
    if (missing(source) || missing(target))
        report(OL$Error, "Source and target images must be given")
    if (!is.nifti(source) || !is.nifti(target))
        report(OL$Error, "Source and target images must be \"nifti\" objects")
    if (!(source@dim_[1] %in% c(2,3,4)))
        report(OL$Error, "Only 2D, 3D or 4D source images may be used")
    if (!(target@dim_[1] %in% c(2,3)))
        report(OL$Error, "Only 2D or 3D target images may be used")
    if (source@dim_[1] == 4 && target@dim_[1] == 2)
        report(OL$Error, "4D to 2D registration cannot be performed")
    if (any(dim(source) < 4) || any(dim(target) < 4))
        report(OL$Error, "Images of fewer than 4 voxels in any dimension cannot be registered")
    if (!is.null(targetMask) && !is.nifti(targetMask))
        report(OL$Error, "Target mask must be NULL or a \"nifti\" object")
    if (any(sapply(list(nLevels,maxIterations,useBlockPercentage,finalInterpolation,verbose), length) != 1))
        report(OL$Error, "Control parameters must all be of unit length")
    if (!(finalInterpolation %in% c(0,1,3)))
        report(OL$Error, "Final interpolation specifier must be 0, 1 or 3")
    
    if (!is.list(initAffine))
        initAffine <- list(initAffine)
    if (!is.null(initAffine[[1]]))
    {
        if (!is.matrix(initAffine[[1]]) || !isTRUE(all.equal(dim(initAffine[[1]]), c(4,4))))
            report(OL$Error, "Specified affine matrix is not valid")
        else if (!is.null(attr(initAffine[[1]],"affineType")) && attr(initAffine[[1]],"affineType") != "niftyreg")
            initAffine <- lapply(initAffine, convertAffine, source=source, target=target, newType="niftyreg")
        
        initAffine <- lapply(initAffine, as.vector, "numeric")
    }
    
    scope <- match.arg(scope)
    
    if (source@dim_[1] == target@dim_[1])
    {
        returnValue <- .Call("reg_aladin_R", .fixTypes(source), .fixTypes(target), scope, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[1]], as.integer(verbose), PACKAGE="RNiftyReg")
        
        dim(returnValue[[1]]) <- dim(target)
        dim(returnValue[[2]]) <- c(4,4)
        attr(returnValue[[2]], "affineType") <- "niftyreg"

        resultImage <- as.nifti(returnValue[[1]], target)
        affine <- list(returnValue[[2]])
        iterations <- list(returnValue[[3]])
    }
    else
    {
        nSourceDims <- source@dim_[1]
        finalDims <- c(dim(target), dim(source)[nSourceDims])
        nReps <- finalDims[length(finalDims)]
        finalArray <- array(0, dim=finalDims)
        affine <- iterations <- vector("list", nReps)
        
        if (length(initAffine) == 1)
            initAffine <- rep(initAffine, nReps)
        else if (length(initAffine) != nReps)
            report(OL$Error, "One initial affine matrix should be provided for each of the ", nReps, " registrations")
        
        for (i in seq_len(nReps))
        {
            if (nSourceDims == 3)
            {
                returnValue <- .Call("reg_aladin_R", .fixTypes(as.nifti(source[,,i],source)), .fixTypes(target), scope, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[i]], as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,i] <- returnValue[[1]]
            }
            else if (nSourceDims == 4)
            {
                returnValue <- .Call("reg_aladin_R", .fixTypes(as.nifti(source[,,,i],source)), .fixTypes(target), scope, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[i]], as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,,i] <- returnValue[[1]]
            }
            
            dim(returnValue[[2]]) <- c(4,4)
            attr(returnValue[[2]], "affineType") <- "niftyreg"
            affine[[i]] <- returnValue[[2]]
            iterations[[i]] <- returnValue[[3]]
        }
        
        resultImage <- as.nifti(finalArray, target)
        resultImage@dim_[nSourceDims+1] <- nReps
    }
    
    resultImage <- .setImageMetadata(resultImage, source, finalInterpolation)
    
    result <- list(image=resultImage, affine=affine, control=NULL, reverseImage=NULL, reverseControl=NULL, iterations=iterations, scope=scope)
    class(result) <- "niftyreg"
    
    return (result)
}

niftyreg.nonlinear <- function (source, target, targetMask = NULL, initAffine = NULL, initControl = NULL, symmetric = FALSE, sourceMask = NULL, nLevels = 3, maxIterations = 300, nBins = 64, bendingEnergyWeight = 0.005, jacobianWeight = 0, inverseConsistencyWeight = 0.01, finalSpacing = c(5,5,5), spacingUnit = c("vox","mm"), finalInterpolation = 3, verbose = FALSE)
{
    if (!require("oro.nifti"))
        report(OL$Error, "The \"oro.nifti\" package is required")
    if (missing(source) || missing(target))
        report(OL$Error, "Source and target images must be given")
    if (!is.nifti(source) || !is.nifti(target))
        report(OL$Error, "Source and target images must be \"nifti\" objects")
    if (!(source@dim_[1] %in% c(2,3,4)))
        report(OL$Error, "Only 2D, 3D or 4D source images may be used")
    if (!(target@dim_[1] %in% c(2,3)))
        report(OL$Error, "Only 2D or 3D target images may be used")
    if (source@dim_[1] == 4 && target@dim_[1] == 2)
        report(OL$Error, "4D to 2D registration cannot be performed")
    if (symmetric && source@dim_[1] != target@dim_[1])
        report(OL$Error, "Source and target images must have the same dimensionality for symmetric registration")
    if (!is.null(targetMask) && !is.nifti(targetMask))
        report(OL$Error, "Target mask must be NULL or a \"nifti\" object")
    if (!is.null(sourceMask) && !is.nifti(sourceMask))
        report(OL$Error, "Source mask must be NULL or a \"nifti\" object")
    if (any(sapply(list(symmetric,nLevels,maxIterations,nBins,bendingEnergyWeight,jacobianWeight,inverseConsistencyWeight,finalInterpolation,verbose), length) != 1))
        report(OL$Error, "Control parameters must all be of unit length")
    if (any(c(bendingEnergyWeight,jacobianWeight,inverseConsistencyWeight) < 0))
        report(OL$Error, "Penalty term weights must be nonnegative")
    if (bendingEnergyWeight + jacobianWeight > 1)
        report(OL$Error, "Penalty term weights cannot add up to more than 1")
    if (!(finalInterpolation %in% c(0,1,3)))
        report(OL$Error, "Final interpolation specifier must be 0, 1 or 3")
    
    # This takes priority over any affine initialisation, if present
    if (!is.list(initControl))
        initControl <- list(initControl)
    if (!is.null(initControl[[1]]))
    {
        if (!is.nifti(initControl[[1]]))
            report(OL$Error, "Initial control point images must be specified as \"nifti\" objects")
        initControl <- lapply(initControl, .fixTypes)
        finalSpacing <- initControl[[1]]@pixdim[2:4]
        spacingUnit <- "mm"
        initAffine <- NULL
    }
    
    if (!is.list(initAffine))
        initAffine <- list(initAffine)
    if (!is.null(initAffine[[1]]))
    {
        if (!is.matrix(initAffine[[1]]) || !isTRUE(all.equal(dim(initAffine[[1]]), c(4,4))))
            report(OL$Error, "Specified affine matrix is not valid")
        else if (!is.null(attr(initAffine[[1]],"affineType")) && attr(initAffine[[1]],"affineType") != "niftyreg")
            initAffine <- lapply(initAffine, convertAffine, source=source, target=target, newType="niftyreg")
        
        initAffine <- lapply(initAffine, as.vector, "numeric")
    }
    
    spacingUnit <- match.arg(spacingUnit)
    if (spacingUnit == "vox")
        finalSpacing <- finalSpacing * abs(target@pixdim[2:4])
    
    if (target@dim_[1] == 2)
    {
        finalSpacing <- c(finalSpacing[1:2], 1)
        controlPointDims <- floor(abs(target@dim_[2:3] * target@pixdim[2:3] / finalSpacing[1:2])) + 5
        controlPointDims <- c(controlPointDims, 1, 1, 2)
        if (symmetric)
        {
            reverseControlPointDims <- floor(abs(source@dim_[2:3] * source@pixdim[2:3] / finalSpacing[1:2])) + 5
            reverseControlPointDims <- c(reverseControlPointDims, 1, 1, 2)
        }
    }
    else
    {
        finalSpacing <- finalSpacing[1:3]
        controlPointDims <- floor(abs(target@dim_[2:4] * target@pixdim[2:4] / finalSpacing)) + 5
        controlPointDims <- c(controlPointDims, 1, 3)
        if (symmetric)
        {
            reverseControlPointDims <- floor(abs(source@dim_[2:4] * source@pixdim[2:4] / finalSpacing)) + 5
            reverseControlPointDims <- c(reverseControlPointDims, 1, 3)
        }
    }
    
    reverseImage <- reverseControl <- NULL
    if (!symmetric)
        sourceMask <- NULL
    
    if (source@dim_[1] == target@dim_[1])
    {
        returnValue <- .Call("reg_f3d_R", .fixTypes(source), .fixTypes(target), as.integer(nLevels), as.integer(maxIterations), as.integer(nBins), as.numeric(bendingEnergyWeight), as.numeric(jacobianWeight), as.numeric(inverseConsistencyWeight), as.numeric(abs(finalSpacing)), as.integer(finalInterpolation), .fixTypes(targetMask), .fixTypes(sourceMask), initAffine[[1]], initControl[[1]], as.integer(symmetric), as.integer(verbose), PACKAGE="RNiftyReg")
        
        dim(returnValue[[1]]) <- dim(target)
        resultImage <- as.nifti(returnValue[[1]], target)
        control <- list(.createControlPointImage(returnValue[[2]], controlPointDims, finalSpacing, returnValue[[3]]))
        iterations <- list(returnValue[[4]])
        
        if (symmetric)
        {
            dim(returnValue[[5]]) <- dim(source)
            reverseImage <- as.nifti(returnValue[[5]], source)
            reverseControl <- list(.createControlPointImage(returnValue[[6]], reverseControlPointDims, finalSpacing, returnValue[[7]]))
        }
    }
    else
    {
        nSourceDims <- source@dim_[1]
        finalDims <- c(dim(target), dim(source)[nSourceDims])
        nReps <- finalDims[length(finalDims)]
        finalArray <- array(0, dim=finalDims)
        control <- iterations <- vector("list", nReps)
        
        if (length(initControl) == 1)
            initControl <- rep(initControl, nReps)
        else if (length(initControl) != nReps)
            report(OL$Error, "One initial control point image should be provided for each of the ", nReps, " registrations")
        
        if (length(initAffine) == 1)
            initAffine <- rep(initAffine, nReps)
        else if (length(initAffine) != nReps)
            report(OL$Error, "One initial affine matrix should be provided for each of the ", nReps, " registrations")

        for (i in seq_len(nReps))
        {
            if (nSourceDims == 3)
            {
                returnValue <- .Call("reg_f3d_R", .fixTypes(as.nifti(source[,,i],source)), .fixTypes(target), as.integer(nLevels), as.integer(maxIterations), as.integer(nBins), as.numeric(bendingEnergyWeight), as.numeric(jacobianWeight), as.numeric(inverseConsistencyWeight), as.numeric(abs(finalSpacing)), as.integer(finalInterpolation), .fixTypes(targetMask), .fixTypes(sourceMask), initAffine[[i]], initControl[[i]], as.integer(symmetric), as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,i] <- returnValue[[1]]
            }
            else if (nSourceDims == 4)
            {
                returnValue <- .Call("reg_f3d_R", .fixTypes(as.nifti(source[,,,i],source)), .fixTypes(target), as.integer(nLevels), as.integer(maxIterations), as.integer(nBins), as.numeric(bendingEnergyWeight), as.numeric(jacobianWeight), as.numeric(inverseConsistencyWeight), as.numeric(abs(finalSpacing)), as.integer(finalInterpolation), .fixTypes(targetMask), .fixTypes(sourceMask), initAffine[[i]], initControl[[i]], as.integer(symmetric), as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,,i] <- returnValue[[1]]
            }
            
            control[[i]] <- .createControlPointImage(returnValue[[2]], controlPointDims, finalSpacing, returnValue[[3]])
            iterations[[i]] <- returnValue[[4]]
        }
        
        resultImage <- as.nifti(finalArray, target)
        resultImage@dim_[nSourceDims+1] <- nReps
    }
    
    resultImage <- .setImageMetadata(resultImage, source, finalInterpolation)
    if (symmetric)
        reverseImage <- .setImageMetadata(reverseImage, target, finalInterpolation)
    
    result <- list(image=resultImage, affine=NULL, control=control, reverseImage=reverseImage, reverseControl=reverseControl, iterations=iterations, scope="nonlinear")
    class(result) <- "niftyreg"
    
    return (result)
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
    
    return (niftyreg.linear(source, target, targetMask=NULL, initAffine=affine, scope="affine", nLevels=0, finalInterpolation=finalInterpolation, verbose=FALSE))
}

applyControlPoints <- function (controlPointImage, source, target, finalInterpolation = 3)
{
    return (niftyreg.nonlinear(source, target, targetMask=NULL, initControl=controlPointImage, symmetric=FALSE, nLevels=0, finalInterpolation=finalInterpolation, verbose=FALSE))
}
