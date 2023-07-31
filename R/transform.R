# For internal use: attach relevant attributes to a transform
# NULL attributes will not be set (as this removes them) unless remove=TRUE is specified
xfmAttrib <- function (transform, source = NULL, target = NULL, ..., remove = FALSE)
{
    if (!is.null(source) && !inherits(source,"niftiImage"))
        source <- asNifti(source, internal=TRUE)
    if (!is.null(target) && !inherits(target,"niftiImage"))
        target <- asNifti(target, internal=TRUE)
    
    attribs <- list(source=source, target=target, ...)
    if (!remove)
        attribs <- attribs[!sapply(attribs,is.null)]
    result <- do.call(structure, c(list(transform),attribs))
    return (result)
}


#' Calculate the deformation field for a transformation
#' 
#' This function is used to calculate the deformation field corresponding to a
#' specified linear or nonlinear transformation. The deformation field gives
#' the location in source image space corresponding to the centre of each voxel
#' in target space. It is used as a common form for linear and nonlinear
#' transformations, and allows them to be visualised.
#' 
#' @param transform A transform, possibly obtained from \code{\link{forward}}
#'   or \code{\link{reverse}}.
#' @param jacobian A logical value: if \code{TRUE}, a Jacobian determinant map
#'   is also calculated and returned in an attribute.
#' @return An \code{"internalImage"} representing the deformation field. If
#'   requested, the Jacobian map is stored in an attribute, which can be
#'   extracted using the \code{\link{jacobian}} accessor function.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg.linear}}, \code{\link{niftyreg.nonlinear}}
#' @export
deformationField <- function (transform, jacobian = TRUE)
{
    if (!isAffine(transform,strict=TRUE) && !isImage(transform,FALSE))
        stop("Specified transformation does not seem to be valid")
    
    return (.Call(C_getDeformationField, transform, isTRUE(jacobian)))
}


#' Extract a Jacobian determinant map
#' 
#' This function extracts the Jacobian determinant map associated with a
#' deformation field.
#' 
#' @param x An R object, probably a deformation field.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{deformationField}}
#' @export
jacobian <- function (x)
{
    return (attr(x, "jacobian"))
}


#' Apply a precomputed transformation
#' 
#' This function allows a precomputed transformation to be applied to a new
#' image or set of points.
#' 
#' Points may be transformed from source to target space exactly under an
#' affine transformation, but nonlinear transformation is inexact. Its accuracy
#' will depend to some extent on the density of the control point grid and the
#' geometry of the deformation in the vicinity of the points of interest.
#' Nevertheless, it should be quite sufficient for most purposes.
#' 
#' The method is to first convert the control points to a deformation field
#' (cf. \code{\link{deformationField}}), which encodes the location of each
#' target space voxel in the source space. The target voxel closest to the
#' requested location is found by searching through this deformation field, and
#' returned if \code{nearest} is \code{TRUE} or it coincides exactly with the
#' requested location. Otherwise, a block of four voxels in each dimension
#' around the point of interest is extracted from the deformation field, and
#' the final location is estimated by local cubic spline regression.
#' 
#' @param transform A transform, possibly obtained from \code{\link{forward}}
#'   or \code{\link{reverse}}.
#' @param x A numeric vector, representing a pixel/voxel location in source
#'   space, or a matrix with rows representing such points, or an image with
#'   the same dimensions as the original source image.
#' @param interpolation A single integer specifying the type of interpolation
#'   to be applied to the final resampled image. May be 0 (nearest neighbour),
#'   1 (trilinear) or 3 (cubic spline). No other values are valid.
#' @param nearest Logical value: if \code{TRUE} and \code{x} contains points,
#'   the nearest voxel centre location in target space will be returned.
#'   Otherwise a more precise subvoxel location will be given.
#' @param internal If \code{FALSE}, the default, the returned image will be
#'   returned as a standard R array. If \code{TRUE}, it will instead be an
#'   object of class \code{"internalImage"}, containing only basic metadata and
#'   a C-level pointer to the full image. (See also \code{\link{readNifti}}.)
#'   This can occasionally be useful to save memory.
#' @return A resampled image or matrix of transformed points.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg.linear}}, \code{\link{niftyreg.nonlinear}}
#' @export
applyTransform <- function (transform, x, interpolation = 3L, nearest = FALSE, internal = FALSE)
{
    source <- attr(transform, "source")
    target <- attr(transform, "target")
    nSourceDim <- ndim(source)
    
    # We only ever return the image, so we don't need the full spectrum of "internal" options
    if (!isTRUE(internal))
        internal <- NA
    
    if (isAffine(transform, strict=TRUE))
    {
        # The argument looks like a suitable image
        if (isImage(x,TRUE) && isTRUE(all.equal(dim(x)[1:nSourceDim],dim(source))))
        {
            result <- niftyreg.linear(x, target, "affine", init=transform, nLevels=0L, interpolation=interpolation, verbose=FALSE, estimateOnly=FALSE, internal=internal)
            return (result$image)
        }
        else if ((is.matrix(x) && ncol(x) == ndim(source)) || length(x) == ndim(source))
        {
            points <- voxelToWorld(x, source, simple=FALSE)
            newPoints <- applyAffine(invertAffine(transform), points)
            newPoints <- worldToVoxel(newPoints, target, simple=FALSE)
            if (nearest)
                newPoints <- round(newPoints)
            return (newPoints)
        }
        else
            stop("Object to transform should be a suitable image or matrix of points")
    }
    else if (isImage(transform, FALSE))
    {
        if (isImage(x,TRUE) && isTRUE(all.equal(dim(x)[1:nSourceDim],dim(source))))
        {
            result <- niftyreg.nonlinear(x, target, init=transform, nLevels=0L, interpolation=interpolation, verbose=FALSE, estimateOnly=FALSE, internal=internal)
            return (result$image)
        }
        else if ((is.matrix(x) && ncol(x) == ndim(source)) || length(x) == ndim(source))
        {
            points <- voxelToWorld(x, source)
            
            if (!is.matrix(points))
                points <- matrix(points, nrow=1)
            
            nDims <- ncol(points)
            if (nDims != ndim(source))
                stop("Dimensionality of points should match the original source image")
            
            result <- .Call(C_transformPoints, transform, points, isTRUE(nearest))
            
            newPoints <- sapply(seq_len(nrow(points)), function(i) {
                if (length(result[[i]]) == nDims)
                    return (result[[i]])
                else
                {
                    data <- as.data.frame(matrix(result[[i]], ncol=2*nDims, byrow=TRUE))
                    if (nDims == 2)
                    {
                        colnames(data) <- c("sx", "sy", "tx", "ty")
                        fit <- lm(cbind(tx,ty) ~ bs(sx) * bs(sy), data=data)
                        return (drop(predict(fit, data.frame(sx=points[i,1],sy=points[i,2]))))
                    }
                    else
                    {
                        colnames(data) <- c("sx", "sy", "sz", "tx", "ty", "tz")
                        fit <- lm(cbind(tx,ty,tz) ~ bs(sx) * bs(sy) * bs(sz), data=data)
                        return (drop(predict(fit, data.frame(sx=points[i,1],sy=points[i,2],sz=points[i,3]))))
                    }
                }
            })
            
            dimnames(newPoints) <- NULL
            newPoints <- drop(t(newPoints))
            return (newPoints)
        }
        else
            stop("Object to transform should be a suitable image or matrix of points")
    }
    else
        stop("Specified transform is not valid")
}


#' Save and load transform objects
#' 
#' These objects save a full transformation object, including source and target
#' image metadata, to a self-contained RDS file, or load it back from such a
#' file.
#' 
#' @param transform A transform, possibly obtained from \code{\link{forward}}
#'   or \code{\link{reverse}}.
#' @param fileName The file name to save to. If \code{NULL}, the serialised
#'   object is returned directly instead.
#' @param x A file name to read from, or a serialised transform object.
#' @return \code{saveTransform} returns a serialised transform object, if no
#'   filename is given; otherwise it is called for its side-effect of writing
#'   to file. \code{loadTransform} returns a deserialised transform object.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{writeAffine}}, \code{\link{readAffine}}
#' @export
saveTransform <- function (transform, fileName = NULL)
{
    source <- niftiHeader(attr(transform, "source"))
    target <- niftiHeader(attr(transform, "target"))
    
    if (isAffine(transform, strict=TRUE))
    {
        transform <- xfmAttrib(transform, remove=TRUE)
        object <- structure(list(transform=transform, source=source, target=target), class="niftyregRDS")
    }
    else if (isImage(transform, FALSE))
    {
        transform <- list(image=as.array(transform), header=niftiHeader(transform), extensions=extensions(transform))
        object <- structure(list(transform=transform, source=source, target=target), class="niftyregRDS")
    }
    else
        stop("Specified transform is not valid")
    
    if (is.null(file))
        return (object)
    else
        saveRDS(object, fileName)
}


#' @rdname saveTransform
#' @export
loadTransform <- function (x)
{
    if (is.character(x))
        object <- readRDS(x)
    else
        object <- x
    
    if (!inherits(object, "niftyregRDS"))
        stop("The specified argument does not contain a serialised transform")
    
    source <- asNifti(object$source, internal=TRUE)
    target <- asNifti(object$target, internal=TRUE)
    
    if (is.list(object$transform))
    {
        transform <- asNifti(object$transform$image, object$transform$header, internal=TRUE)
        extensions(transform) <- object$transform$extensions
        return (xfmAttrib(transform, source, target))
    }
    else
        return (asAffine(object$transform, source, target))
}


#' Apply simple transformations
#' 
#' These functions allow simple transformations to be applied quickly, or in a
#' chosen order. They represent simplified interfaces to the
#' \code{\link{buildAffine}} and \code{\link{applyTransform}} functions, and
#' are compatible with the chaining operator from the popular \code{magrittr}
#' package (although performing one single transformation may be preferable).
#' 
#' @inheritParams buildAffine
#' @param source A 2D or 3D image, in the sense of \code{\link{isImage}}.
#' @param ... Additional arguments to \code{\link{applyTransform}}.
#' @return The transformed image.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{buildAffine}}, \code{\link{applyTransform}}
#' @export
translate <- function (source, translation, ...)
{
    xfm <- buildAffine(translation=translation, source=source)
    applyTransform(xfm, source, ...)
}


#' @rdname translate
#' @export
rescale <- function (source, scales, anchor = c("none","origin","centre","center"), ...)
{
    anchor <- match.arg(anchor)
    xfm <- buildAffine(scales=scales, source=source, anchor=anchor)
    applyTransform(xfm, source, ...)
}


#' @rdname translate
#' @export
skew <- function (source, skews, anchor = c("none","origin","centre","center"), ...)
{
    anchor <- match.arg(anchor)
    xfm <- buildAffine(skews=skews, source=source, anchor=anchor)
    applyTransform(xfm, source, ...)
}


#' @rdname translate
#' @export
rotate <- function (source, angles, anchor = c("none","origin","centre","center"), ...)
{
    anchor <- match.arg(anchor)
    xfm <- buildAffine(angles=angles, source=source, anchor=anchor)
    applyTransform(xfm, source, ...)
}


#' Calculate a half transformation
#' 
#' This function calculates the half-way transformation corresponding to its
#' argument. Applying this transformation results in points or images in a
#' space halfway between the original source and target images, which can be a
#' useful common space in some applications.
#' 
#' @param transform A transform, possibly obtained from \code{\link{forward}}
#'   or \code{\link{reverse}}.
#' @return The half-way transform, in a similar format to \code{transform}.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg.linear}}, \code{\link{niftyreg.nonlinear}}
#' @export
halfTransform <- function (transform)
{
    invisible (.Call(C_halfTransform, transform))
}


#' Compose transformations
#' 
#' Compute the composition of two or more transforms, the single transform that
#' combines their effects in order.
#' 
#' @param ... Affine or nonlinear transforms, possibly obtained from
#'   \code{\link{forward}} or \code{\link{reverse}}.
#' @return The composed transform. If all arguments are affines then the result
#'   will also be an affine; otherwise it will be a deformation field.
#' 
#' @note The source image for the composed transform is generally the source
#'   image from the first transform, and the target is the target image from
#'   the second transform. However, the target image attached to half
#'   transforms (as calculated by \code{\link{halfTransform}}) generally has a
#'   modified xform, compared to the original target. Therefore, composing a
#'   half transform with itself may not be exactly equivalent to the original.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg.linear}}, \code{\link{niftyreg.nonlinear}},
#'   \code{\link{deformationField}}
#' @export
composeTransforms <- function (...)
{
    composePair <- function(t1,t2) .Call(C_composeTransforms, t1, t2)
    invisible (Reduce(composePair, list(...)))
}
