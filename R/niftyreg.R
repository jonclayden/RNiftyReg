niftyreg <- function (source, target, targetMask = NULL, initAffine = NULL, scope = c("affine","rigid"), nLevels = 3, maxIterations = 5, useBlockPercentage = 50, verbose = FALSE)
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
    if (any(sapply(list(nLevels,maxIterations,useBlockPercentage,verbose), length) != 1))
        stop("Control parameters must all be of unit length")
    
    if (!is.null(initAffine))
    {
        if (!is.matrix(initAffine) || !isTRUE(all.equal(dim(initAffine), c(4,4))))
            stop("Specified affine matrix is not valid")
        else if (!is.null(attr(initAffine,"affineType")) && attr(initAffine,"affineType") == "fsl")
            initAffine <- convertAffine(initAffine, source, target, newType="niftyreg")
        
        initAffine <- as.vector(initAffine, "numeric")
    }
    
    scope <- match.arg(scope)
    
    returnValue <- .Call("reg_aladin", source, target, scope, as.integer(nLevels), as.integer(maxIterations), as.integer(useBlockPercentage), targetMask, initAffine, as.integer(verbose), PACKAGE="RNiftyReg")
    
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
