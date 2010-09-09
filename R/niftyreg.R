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
            stop("RNiftyReg does not support the \"", datatypeName, "\" image data type")
        else if (datatypeName %in% doubleDatatypes)
            storage.mode(image@.Data) <- "double"
        else
            storage.mode(image@.Data) <- "integer"
        
        return (image)
    }
}

niftyreg <- function (source, target, targetMask = NULL, initAffine = NULL, scope = c("affine","rigid"), nLevels = 3, maxIterations = 5, useBlockPercentage = 50, finalInterpolation = 3, verbose = FALSE)
{
    if (!require("oro.nifti"))
        stop("The \"oro.nifti\" package is required")
    if (missing(source) || missing(target))
        stop("Source and target images must be given")
    if (!is.nifti(source) || !is.nifti(target))
        stop("Source and target images must be \"nifti\" objects")
    if (source@dim_[1] != 3 || target@dim_[1] != 3)
        stop("Only 3D source and target images may be used at present")
    if (!is.null(targetMask) && !is.nifti(targetMask))
        stop("Target mask must be NULL or a \"nifti\" object")
    if (any(sapply(list(nLevels,maxIterations,useBlockPercentage,finalInterpolation,verbose), length) != 1))
        stop("Control parameters must all be of unit length")
    if (!(finalInterpolation %in% c(0,1,3)))
        stop("Final interpolation specifier must be 0, 1 or 3")
    
    if (!is.null(initAffine))
    {
        if (!is.matrix(initAffine) || !isTRUE(all.equal(dim(initAffine), c(4,4))))
            stop("Specified affine matrix is not valid")
        else if (!is.null(attr(initAffine,"affineType")) && attr(initAffine,"affineType") == "fsl")
            initAffine <- convertAffine(initAffine, source, target, newType="niftyreg")
        
        initAffine <- as.vector(initAffine, "numeric")
    }
    
    scope <- match.arg(scope)
    
    returnValue <- .Call("reg_aladin", .fixTypes(source), .fixTypes(target), scope, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), as.integer(finalInterpolation), .fixTypes(targetMask), initAffine, as.integer(verbose), PACKAGE="RNiftyReg")
    
    dim(returnValue[[1]]) <- dim(target)
    dim(returnValue[[2]]) <- c(4,4)
    attr(returnValue[[2]], "affineType") <- "niftyreg"
    
    resultImage <- as.nifti(returnValue[[1]], target)
    resultImage@cal_min <- min(resultImage@.Data)
    resultImage@cal_max <- max(resultImage@.Data)
    resultImage@scl_slope <- source@scl_slope
    resultImage@scl_inter <- source@scl_inter
    resultImage@datatype <- source@datatype
    resultImage@bitpix <- as.numeric(source@bitpix)
    
    result <- list(image=resultImage, affine=returnValue[[2]], scope=scope)
    class(result) <- "niftyreg"
    
    return (result)
}
