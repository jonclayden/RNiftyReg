isImage <- function (object, unsure = NA)
{
    if (any(c("nifti","internalImage") %in% class(object)))
        return (TRUE)
    else if (!is.null(attr(object, ".nifti_image_ptr")))
        return (TRUE)
    else if (is.null(dim(object)))
        return (FALSE)
    else if (isAffine(object))
        return (FALSE)
    else
        return (unsure)
}

dim.internalImage <- function (x)
{
    return (attr(x, "imagedim"))
}

as.array.internalImage <- function (x, ...)
{
    return (.Call("pointerToArray", x, PACKAGE="RNiftyReg"))
}

print.internalImage <- function (x, ...)
{
    dim <- attr(x, "imagedim")
    ndim <- length(dim)
    pixdim <- attr(x, "pixdim")
    pixunits <- attr(x, "pixunits")
    cat(paste0("Internal image: \"", x, "\"\n"))
    cat(paste("-", paste(dim,collapse=" x "), ifelse(ndim>2,"voxels\n","pixels\n")))
    
    if (!is.null(pixdim))
    {
        spaceUnit <- grep("m$", pixunits, perl=TRUE, value=TRUE)
        cat(paste("-", paste(signif(pixdim[1:min(3,ndim)],4),collapse=" x ")))
        if (length(spaceUnit) > 0)
            cat(paste0(" ", spaceUnit[1]))
        
        if (ndim > 3)
        {
            timeUnit <- grep("s$", pixunits, perl=TRUE, value=TRUE)
            cat(paste(" x", signif(pixdim[4],4)))
            if (length(timeUnit) > 0)
                cat(paste0(" ", timeUnit[1]))
        }
        if (ndim > 4)
            cat(paste(" x", paste(signif(pixdim[5:ndim],4),collapse=" x ")))
        
        cat(paste(" per", ifelse(ndim>2,"voxel\n","pixel\n")))
    }
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

pixunits <- function (object)
{
    if (!is.null(attr(object, "pixunits")))
        return (attr(object, "pixunits"))
    else
        return ("Unknown")
}
