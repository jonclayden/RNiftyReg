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
    
    return (.Call("getDeformationField", transform, isTRUE(jacobian), PACKAGE="RNiftyReg"))
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
            
            result <- .Call("transformPoints", transform, points, isTRUE(nearest), PACKAGE="RNiftyReg")
            
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
    invisible (.Call("halfTransform", transform, PACKAGE="RNiftyReg"))
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
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg.linear}}, \code{\link{niftyreg.nonlinear}},
#'   \code{\link{deformationField}}
#' @export
composeTransforms <- function (...)
{
    composePair <- function(t1,t2) .Call("composeTransforms", t1, t2, PACKAGE="RNiftyReg")
    invisible (Reduce(composePair, list(...)))
}
