#' Read a NIfTI-1 format file
#' 
#' This function reads one or more NIfTI-1 files into R, using the standard
#' NIfTI-1 C library.
#' 
#' @param file A character vector of file names.
#' @param internal Logical value. If \code{FALSE} (the default), an array
#'   containing the image pixel or voxel values will be returned. If
#'   \code{TRUE}, the return value will be an object of class
#'   \code{"internalImage"}, which contains only minimal metadata about the
#'   image. Either way, the return value has an attribute which points to a
#'   C data structure containing the full image.
#' @return An array or \code{"internalImage"} object.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @references The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}).
#' @export
readNifti <- function (file, internal = FALSE)
{
    if (!is.character(file))
        stop("File name(s) must be specified in a character vector")
    if (length(file) == 0)
        stop("File name vector is empty")
    else if (length(file) > 1)
        lapply(file, function(f) .Call("readNifti", f, internal, PACKAGE="RNiftyReg"))
    else
        .Call("readNifti", file, internal, PACKAGE="RNiftyReg")
}


#' Write a NIfTI-1 format file
#' 
#' This function writes an image to NIfTI-1 format, using the standard NIfTI-1
#' C library.
#' 
#' @param image An image, in any acceptable form (see \code{\link{isImage}}).
#' @param file A character string containing a file name.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @references The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}).
#' @export
writeNifti <- function (image, file)
{
    .Call("writeNifti", image, file, PACKAGE="RNiftyReg")
}
