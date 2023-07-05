#' Read a NIfTI format file
#' 
#' This function reads files in NIfTI-1 or NIfTI-2 format into R, using the
#' standard NIfTI C library. It extends the equivalent function from the
#' \code{RNifti} package with source and target image parameters.
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
#' @export
readNifti <- function (file, source = NULL, target = NULL, internal = FALSE)
{
    image <- RNifti::readNifti(file, internal=internal)
    return (xfmAttrib(image, source, target))
}
