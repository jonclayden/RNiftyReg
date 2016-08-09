#' Test whether an object represents an image
#' 
#' This function tried to determine whether an object is an image that the
#' package knows how to handle. If its class is \code{"nifti"},
#' \code{"niftiImage"}, \code{"internalImage"} or \code{"MriImage"}, then the
#' result is always \code{TRUE}. Likewise if it has an internal image pointer.
#' If it has no \code{dim} attribute, or looks like an affine matrix, then the
#' result is \code{FALSE}. Otherwise the value of the \code{unsure} argument
#' is returned.
#' 
#' @param object An R object.
#' @param unsure The value to return if the function can't tell whether or not
#'   the \code{object} is an image.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
isImage <- function (object, unsure = NA)
{
    if (any(c("nifti","niftiImage","internalImage","MriImage") %in% class(object)))
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
