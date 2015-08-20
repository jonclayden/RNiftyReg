#' Read a NIfTI-1 format file
#' 
#' This function reads one or more NIfTI-1 files into R, using the standard
#' NIfTI-1 C library.
#' 
#' @param file A character vector of file names.
#' @param source,target If the specified \code{file} contains a transformation,
#'   these parameters can be used to specify the associated source and target
#'   images, which are stored in attributes of the same name. Only used if
#'   \code{file} is of unit length.
#' @param internal Logical value. If \code{FALSE} (the default), an array
#'   of class \code{"niftiImage"}, containing the image pixel or voxel values,
#'   will be returned. If \code{TRUE}, the return value will be an object of
#'   class \code{"internalImage"}, which contains only minimal metadata about
#'   the image. Either way, the return value has an attribute which points to a
#'   C data structure containing the full image.
#' @return An array or internal image, with class \code{"niftiImage"}, and
#'   possibly also \code{"internalImage"}.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{writeNifti}}
#' @references The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}).
#' @export
readNifti <- function (file, source = NULL, target = NULL, internal = FALSE)
{
    if (!is.character(file))
        stop("File name(s) must be specified in a character vector")
    if (length(file) == 0)
        stop("File name vector is empty")
    else if (length(file) > 1)
        lapply(file, function(f) .Call("readNifti", path.expand(f), internal, PACKAGE="RNiftyReg"))
    else
    {
        image <- .Call("readNifti", path.expand(file), internal, PACKAGE="RNiftyReg")
        if (!is.null(source))
            attr(image, "source") <- .Call("retrieveImage", source, PACKAGE="RNiftyReg")
        if (!is.null(target))
            attr(image, "target") <- .Call("retrieveImage", target, PACKAGE="RNiftyReg")
        return (image)
    }
}


#' Write a NIfTI-1 format file
#' 
#' This function writes an image to NIfTI-1 format, using the standard NIfTI-1
#' C library.
#' 
#' @param image An image, in any acceptable form (see \code{\link{isImage}}).
#' @param file A character string containing a file name.
#' @param template An optional template object to derive NIfTI-1 properties
#'   from. Passed to \code{\link{updateNifti}} if \code{image} is an array.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{readNifti}}, \code{\link{updateNifti}}
#' @references The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}).
#' @export
writeNifti <- function (image, file, template = NULL)
{
    if (is.array(image) && !is.null(template))
        image <- updateNifti(image, template)
    
    invisible(.Call("writeNifti", image, file, PACKAGE="RNiftyReg"))
}


#' Update an internal NIfTI-1 object using a template
#' 
#' This function adds or updates the internal NIfTI-1 object for an array,
#' using metadata from the template. The dimensions and, if available, pixel
#' dimensions, from the \code{image} will replace those from the template.
#' 
#' @param image A numeric array.
#' @param template An image, in any acceptable form (see
#'   \code{\link{isImage}}). The default of \code{NULL} will have no effect.
#' @return A copy of the original \code{image}, with its internal image
#'   attribute set or updated appropriately.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @export
updateNifti <- function (image, template = NULL)
{
    .Call("updateNifti", image, template, PACKAGE="RNiftyReg")
}


#' Dump the contents of an internal NIfTI-1 object
#' 
#' This function extracts the contents of an internal NIfTI-1 object into an R
#' list. No processing is done to the elements.
#' 
#' @param image An image, in any acceptable form (see \code{\link{isImage}}).
#' @return A list with named components corresponding to the elements in a raw
#'   NIfTI-1 file.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @references The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}).
#' @export
dumpNifti <- function (image)
{
    .Call("dumpNifti", image, PACKAGE="RNiftyReg")
}
