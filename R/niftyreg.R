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

niftyreg.linear <- function (source, target, targetMask = NULL, initAffine = NULL, scope = c("affine","rigid"), nLevels = 3, maxIterations = 5, useBlockPercentage = 50, finalInterpolation = 3, verbose = FALSE, interpolationPrecision = NULL)
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
    if (length(dim(source)) - length(dim(target)) > 1)
        report(OL$Error, "The source image may not have more than one extra dimension")
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
    
    if (!is.null(interpolationPrecision))
        interpolationPrecision <- match.arg(interpolationPrecision, c("source","single","double"))
    else if (finalInterpolation == 0)
        interpolationPrecision <- "source"
    else
        interpolationPrecision <- "single"
    
    if (source@dim_[1] == target@dim_[1])
    {
        returnValue <- .Call("reg_aladin_R", .fixTypes(source), .fixTypes(target), scope, interpolationPrecision, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[1]], as.integer(verbose), PACKAGE="RNiftyReg")
        
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
                returnValue <- .Call("reg_aladin_R", .fixTypes(as.nifti(source[,,i],source)), .fixTypes(target), scope, interpolationPrecision, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[i]], as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,i] <- returnValue[[1]]
            }
            else if (nSourceDims == 4)
            {
                returnValue <- .Call("reg_aladin_R", .fixTypes(as.nifti(source[,,,i],source)), .fixTypes(target), scope, interpolationPrecision, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[i]], as.integer(verbose), PACKAGE="RNiftyReg")
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
    
    resultImage@cal_min <- min(resultImage@.Data)
    resultImage@cal_max <- max(resultImage@.Data)
    resultImage@scl_slope <- source@scl_slope
    resultImage@scl_inter <- source@scl_inter
    
    resultImage@datatype <- switch(interpolationPrecision, source=source@datatype, single=16L, double=64L)
    resultImage@data_type <- convert.datatype(resultImage@datatype)
    resultImage@bitpix <- switch(interpolationPrecision, source=as.numeric(source@bitpix), single=32, double=64)
    
    result <- list(image=resultImage, affine=affine, control=NULL, iterations=iterations, scope=scope)
    class(result) <- "niftyreg"
    
    return (result)
}

niftyreg.nonlinear <- function (source, target, targetMask = NULL, initAffine = NULL, initControl = NULL, nLevels = 3, maxIterations = 300, nBins = 64, bendingEnergyWeight = 0.01, jacobianWeight = 0, finalSpacing = c(5,5,5), spacingUnit = c("vox","mm"), finalInterpolation = 3, verbose = FALSE, interpolationPrecision = NULL)
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
    if (length(dim(source)) - length(dim(target)) > 1)
        report(OL$Error, "The source image may not have more than one extra dimension")
    if (!is.null(targetMask) && !is.nifti(targetMask))
        report(OL$Error, "Target mask must be NULL or a \"nifti\" object")
    if (any(sapply(list(nLevels,maxIterations,nBins,bendingEnergyWeight,jacobianWeight,finalInterpolation,verbose), length) != 1))
        report(OL$Error, "Control parameters must all be of unit length")
    if (any(c(bendingEnergyWeight,jacobianWeight) < 0))
        report(OL$Error, "Penalty term weights must be nonnegative")
    if (bendingEnergyWeight + jacobianWeight > 1)
        report(OL$Error, "Penalty term weights cannot add up to more than 1")
    if (length(finalSpacing) != 3)
        report(OL$Error, "Final spacing must be specified as a numeric 3-vector")
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
        finalSpacing <- finalSpacing * target@pixdim[2:4]
    controlPointDims <- floor(abs(target@dim_[2:4] * target@pixdim[2:4] / finalSpacing)) + 5
    
    if (!is.null(interpolationPrecision))
        interpolationPrecision <- match.arg(interpolationPrecision, c("source","single","double"))
    else if (finalInterpolation == 0)
        interpolationPrecision <- "source"
    else
        interpolationPrecision <- "single"
    
    if (source@dim_[1] == target@dim_[1])
    {
        returnValue <- .Call("reg_f3d_R", .fixTypes(source), .fixTypes(target), interpolationPrecision, as.integer(nLevels), as.integer(maxIterations), as.integer(nBins), as.numeric(bendingEnergyWeight), as.numeric(jacobianWeight), as.numeric(abs(finalSpacing)), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[1]], initControl[[1]], as.integer(verbose), PACKAGE="RNiftyReg")
        
        dim(returnValue[[1]]) <- dim(target)
        dim(returnValue[[2]]) <- c(controlPointDims,1,3)
        
        xform <- returnValue[[4]]

        resultImage <- as.nifti(returnValue[[1]], target)
        control <- list(new("nifti", .Data=returnValue[[2]], dim_=c(5,controlPointDims,1,3,1,1), datatype=64L, bitpix=64, pixdim=c(xform[9],finalSpacing,1,1,0,0), xyzt_units=0, qform_code=xform[1], sform_code=xform[2], quatern_b=xform[3], quatern_c=xform[4], quatern_d=xform[5], qoffset_x=xform[6], qoffset_y=xform[7], qoffset_z=xform[8], srow_x=xform[10:13], srow_y=xform[14:17], srow_z=xform[18:21], cal_min=min(returnValue[[2]]), cal_max=max(returnValue[[2]])))
        iterations <- list(returnValue[[3]])
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
                returnValue <- .Call("reg_f3d_R", .fixTypes(as.nifti(source[,,i],source)), .fixTypes(target), interpolationPrecision, as.integer(nLevels), as.integer(maxIterations), as.integer(nBins), as.numeric(bendingEnergyWeight), as.numeric(jacobianWeight), as.numeric(abs(finalSpacing)), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[i]], initControl[[i]], as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,i] <- returnValue[[1]]
            }
            else if (nSourceDims == 4)
            {
                returnValue <- .Call("reg_f3d_R", .fixTypes(as.nifti(source[,,,i],source)), .fixTypes(target), interpolationPrecision, as.integer(nLevels), as.integer(maxIterations), as.integer(nBins), as.numeric(bendingEnergyWeight), as.numeric(jacobianWeight), as.numeric(abs(finalSpacing)), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine[[i]], initControl[[i]], as.integer(verbose), PACKAGE="RNiftyReg")
                finalArray[,,,i] <- returnValue[[1]]
            }
            
            dim(returnValue[[2]]) <- c(controlPointDims,1,3)

            xform <- returnValue[[4]]

            control[[i]] <- new("nifti", .Data=returnValue[[2]], dim_=c(5,controlPointDims,1,3,1,1), datatype=64L, bitpix=64, pixdim=c(xform[9],finalSpacing,1,1,0,0), xyzt_units=0, qform_code=xform[1], sform_code=xform[2], quatern_b=xform[3], quatern_c=xform[4], quatern_d=xform[5], qoffset_x=xform[6], qoffset_y=xform[7], qoffset_z=xform[8], srow_x=xform[10:13], srow_y=xform[14:17], srow_z=xform[18:21], cal_min=min(returnValue[[2]]), cal_max=max(returnValue[[2]]))
            iterations[[i]] <- returnValue[[3]]
        }
        
        resultImage <- as.nifti(finalArray, target)
        resultImage@dim_[nSourceDims+1] <- nReps
    }
    
    resultImage@cal_min <- min(resultImage@.Data)
    resultImage@cal_max <- max(resultImage@.Data)
    resultImage@scl_slope <- source@scl_slope
    resultImage@scl_inter <- source@scl_inter
    
    resultImage@datatype <- switch(interpolationPrecision, source=source@datatype, single=16L, double=64L)
    resultImage@data_type <- convert.datatype(resultImage@datatype)
    resultImage@bitpix <- switch(interpolationPrecision, source=as.numeric(source@bitpix), single=32, double=64)
    
    result <- list(image=resultImage, affine=NULL, control=control, iterations=iterations, scope="nonlinear")
    class(result) <- "niftyreg"
    
    return (result)
}

applyAffine <- function (affine, source, target, affineType = NULL, finalInterpolation = 3, interpolationPrecision = NULL)
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
    
    return (niftyreg.linear(source, target, targetMask=NULL, initAffine=affine, scope="affine", nLevels=0, finalInterpolation=finalInterpolation, interpolationPrecision=interpolationPrecision, verbose=FALSE))
}

applyControlPoints <- function (controlPointImage, source, target, finalInterpolation = 3, interpolationPrecision = NULL)
{
    return (niftyreg.nonlinear(source, target, targetMask=NULL, initControl=controlPointImage, nLevels=0, finalInterpolation=finalInterpolation, interpolationPrecision=interpolationPrecision, verbose=FALSE))
}
