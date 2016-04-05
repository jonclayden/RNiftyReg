#' Two and three dimensional image registration
#' 
#' The \code{niftyreg} function performs linear or nonlinear registration for
#' two and three dimensional images. 4D images may also be registered
#' volumewise to a 3D image, or 3D images slicewise to a 2D image. This
#' function is a common wrapper for \code{\link{niftyreg.linear}} and
#' \code{\link{niftyreg.nonlinear}}.
#' 
#' @param source The source image, an object of class \code{"nifti"} or
#'   \code{"internalImage"}, or a plain array, or a NIfTI-1 filename. Must have
#'   2, 3 or 4 dimensions.
#' @param target The target image, an object of class \code{"nifti"} or
#'   \code{"internalImage"}, or a plain array, or a NIfTI-1 filename. Must have
#'   2 or 3 dimensions.
#' @param scope A string describing the scope, or number of degrees of freedom
#'   (DOF), of the registration. The currently supported values are
#'   \code{"affine"} (12 DOF), \code{"rigid"} (6 DOF) or \code{"nonlinear"}
#'   (high DOF, with the exact number depending on the image sizes).
#' @param init Transformation(s) to be used for initialisation, which may be
#'   \code{NULL}, for no initialisation, or an affine matrix or control point
#'   image (nonlinear only). For multiple registration, where the source image
#'   has one more dimension than the target, this may also be a list whose
#'   components are likewise \code{NULL} or a suitable initial transform.
#' @param sourceMask An optional mask image in source space, whose nonzero
#'   region will be taken as the region of interest for the registration.
#'   Ignored when \code{symmetric} is \code{FALSE}.
#' @param targetMask An optional mask image in target space, whose nonzero
#'   region will be taken as the region of interest for the registration.
#' @param symmetric Logical value. Should forward and reverse transformations
#'   be estimated simultaneously?
#' @param interpolation A single integer specifying the type of interpolation
#'   to be applied to the final resampled image. May be 0 (nearest neighbour),
#'   1 (trilinear) or 3 (cubic spline). No other values are valid.
#' @param estimateOnly Logical value: if \code{TRUE}, transformations will be
#'   estimated, but images will not be resampled.
#' @param sequentialInit If \code{TRUE} and \code{source} has higher
#'   dimensionality than \code{target}, transformations which are not
#'   explicitly initialised will begin from the result of the previous
#'   registration.
#' @param internal If \code{FALSE}, the default, the returned image will be
#'   returned as a standard R array. If \code{TRUE}, it will instead be an
#'   object of class \code{"internalImage"}, containing only basic metadata and
#'   a C-level pointer to the full image. (See also \code{\link{readNifti}}.)
#'   This can occasionally be useful to save memory.
#' @param ... Further arguments to \code{\link{niftyreg.linear}} or
#'   \code{\link{niftyreg.nonlinear}}.
#' @param x A \code{"niftyreg"} object.
#' @return A list of class \code{"niftyreg"} with components:
#'   \describe{
#'     \item{image}{An array or internal image representing the registered and
#'       resampled \code{source} image in the space of the \code{target} image.
#'       This element is \code{NULL} if the \code{estimateOnly} parameter is
#'       \code{TRUE}.}
#'     \item{forwardTransforms}{A list of (linear or nonlinear) transformations
#'       from source to target space.}
#'     \item{reverseTransforms}{A list of (linear or nonlinear) transformations
#'       from target to source space.}
#'     \item{iterations}{A list of integer vectors, giving the number of
#'       iterations completed at each ``level'' of the algorithm. Note that for
#'       the first level of the linear algorithm specifically, twice the
#'       specified number of iterations is allowed.}
#'     \item{source}{An internal representation of the source image for each
#'       registration.}
#'     \item{target}{An internal representation of the target image.}
#'   }
#'   The \code{as.array} method for this class returns the \code{image}
#'   element.
#' 
#' @note If substantial parts of the target image are zero-valued, for example
#'   because the target image has been brain-extracted, it can be useful to
#'   pass it as a target mask as well as the target image, viz.
#'   \code{niftyreg(source, target, targetMask=target)}.
#' 
#' @examples
#' \dontrun{
#' source <- readNifti(system.file("extdata", "epi_t2.nii.gz",
#'   package="RNiftyReg"))
#' target <- readNifti(system.file("extdata", "flash_t1.nii.gz",
#'   package="RNiftyReg"))
#' 
#' result <- niftyreg(source, target, scope="affine")
#' }
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg.linear}} and \code{\link{niftyreg.nonlinear}},
#'   which do most of the work. Also, \code{\link{forward}} and
#'   \code{\link{reverse}} to extract transformations, and
#'   \code{\link{applyTransform}} to apply them to new images or points.
#' @references Please see \code{\link{niftyreg.linear}} or
#' \code{\link{niftyreg.nonlinear}} for references relating to each type of
#' registration.
#' @export
niftyreg <- function (source, target, scope = c("affine","rigid","nonlinear"), init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, interpolation = 3L, estimateOnly = FALSE, sequentialInit = FALSE, internal = FALSE, ...)
{
    if (missing(source) || missing(target))
        stop("Source and target images must be given")
    
    scope <- match.arg(scope)
    if (scope == "nonlinear")
        niftyreg.nonlinear(source, target, init, sourceMask, targetMask, interpolation=interpolation, symmetric=symmetric, estimateOnly=estimateOnly, sequentialInit=sequentialInit, internal=internal, ...)
    else
        niftyreg.linear(source, target, scope, init, sourceMask, targetMask, interpolation=interpolation, symmetric=symmetric, estimateOnly=estimateOnly, sequentialInit=sequentialInit, internal=internal, ...)
}


#' Two and three dimensional linear image registration
#' 
#' The \code{niftyreg.linear} function performs linear registration for two and
#' three dimensional images. 4D images may also be registered volumewise to a
#' 3D image, or 3D images slicewise to a 2D image. Rigid-body (6 degrees of
#' freedom) and affine (12 degrees of freedom) registration can currently be
#' performed.
#' 
#' This function performs the dual operations of finding a transformation to
#' optimise image alignment, and resampling the source image into the space of
#' the target image.
#' 
#' The algorithm is based on a block-matching approach and Least Trimmed
#' Squares (LTS) fitting. Firstly, the block matching provides a set of
#' corresponding points between a target and a source image. Secondly, using
#' this set of corresponding points, the best rigid or affine transformation is
#' evaluated. This two-step loop is repeated until convergence to the best
#' transformation is achieved.
#' 
#' In the NiftyReg implementation, normalised cross-correlation between the
#' target and source blocks is used to evaluate correspondence. The block width
#' is constant and has been set to 4 voxels. A coarse-to-fine approach is used,
#' where the registration is first performed on down-sampled images (using a
#' Gaussian filter to resample images), and finally performed on full
#' resolution images.
#' 
#' The source image may have 2, 3 or 4 dimensions, and the target 2 or 3. The
#' dimensionality of the target image determines whether 2D or 3D registration
#' is applied, and source images with one more dimension than the target (i.e.
#' 4D to 3D, or 3D to 2D) will be registered volumewise or slicewise, as
#' appropriate. In the latter case the last dimension of the resulting image is
#' taken from the source image, while all other dimensions come from the
#' target. One affine matrix is returned for each registration performed.
#'
#' @inheritParams niftyreg 
#' @param nLevels A single integer specifying the number of levels of the
#'   algorithm that should be applied. If zero, no optimisation will be
#'   performed, and the final affine matrix will be the same as its
#'   initialisation value.
#' @param maxIterations A single integer specifying the maximum number of
#'   iterations to be used within each level. Fewer iterations may be used if a
#'   convergence test deems the process to have completed.
#' @param useBlockPercentage A single integer giving the percentage of blocks
#'   to use for calculating correspondence at each step of the algorithm. The
#'   blocks with the highest intensity variance will be chosen.
#' @param verbose A single logical value: if \code{TRUE}, the code will give
#'   some feedback on its progress; otherwise, nothing will be output while the
#'   algorithm runs. Run time can be seconds or more, depending on the size and
#'   dimensionality of the images.
#' @return See \code{\link{niftyreg}}.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg}}, which can be used as an interface to this
#'   function, and \code{\link{niftyreg.nonlinear}} for nonlinear registration.
#'   Also, \code{\link{forward}} and \code{\link{reverse}} to extract
#'   transformations, and \code{\link{applyTransform}} to apply them to new
#'   images or points.
#' @references The algorithm used by this function is described in the
#' following publication.
#' 
#' M. Modat, D.M. Cash, P. Daga, G.P. Winston, J.S. Duncan & S. Ourselin
#' (2014). Global image registration using a symmetric block-matching approach.
#' Journal of Medical Imaging 1(2):024003.
#' @export
niftyreg.linear <- function (source, target, scope = c("affine","rigid"), init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, nLevels = 3L, maxIterations = 5L, useBlockPercentage = 50L, interpolation = 3L, verbose = FALSE, estimateOnly = FALSE, sequentialInit = FALSE, internal = FALSE)
{
    if (missing(source) || missing(target))
        stop("Source and target images must be given")
    
    source <- .Call("retrieveImage", source, PACKAGE="RNiftyReg")
    target <- .Call("retrieveImage", target, PACKAGE="RNiftyReg")
    nSourceDim <- ndim(source)
    nTargetDim <- ndim(target)
    
    if (!(interpolation %in% c(0,1,3)))
        stop("Final interpolation specifier must be 0, 1 or 3")
    
    scope <- match.arg(scope)
    nReps <- ifelse(nSourceDim > nTargetDim, dim(source)[nSourceDim], 1L)
    
    if (!is.list(init))
        init <- list(init)
    if (length(init) != nReps)
    {
        if (sequentialInit)
            init <- c(init, rep(list(NULL),nReps-length(init)))
        else
            init <- rep(init, length.out=nReps)
    }
    init <- lapply(init, function(x) {
        if (!is.null(x) && !isAffine(x))
            stop("Linear registration can only be initialised with an affine matrix")
        else
            return (x)
    })
    
    result <- .Call("regLinear", source, target, ifelse(scope=="affine",1L,0L), symmetric, nLevels, maxIterations, useBlockPercentage, interpolation, sourceMask, targetMask, init, verbose, estimateOnly, sequentialInit, internal, PACKAGE="RNiftyReg")
    class(result) <- "niftyreg"
    
    return (result)
}


#' Two and three dimensional nonlinear image registration
#' 
#' The \code{niftyreg.nonlinear} function performs nonlinear registration for
#' two and three dimensional images. 4D images may also be registered
#' volumewise to a 3D image, or 3D images slicewise to a 2D image. The warping
#' is based on free-form deformations, parameterised using an image of control
#' points.
#' 
#' This function performs the dual operations of finding a transformation to
#' optimise image alignment, and resampling the source image into the space of
#' the target image (and vice-versa, if \code{symmetric} is \code{TRUE}).
#' Unlike \code{\link{niftyreg.linear}}, this transformation is nonlinear, and
#' the degree of deformation may vary across the image.
#' 
#' The nonlinear warping is based on free-form deformations. A lattice of
#' equally-spaced control points is defined over the target image, each of
#' which can be moved to locally modify the mapping to the source image. In
#' order to assess the quality of the warping between the two images, an
#' objective function based on the normalised mutual information is used, with
#' penalty terms based on the bending energy or the squared log of the Jacobian
#' determinant. The objective function value is optimised using a conjugate
#' gradient scheme.
#' 
#' The source image may have 2, 3 or 4 dimensions, and the target 2 or 3. The
#' dimensionality of the target image determines whether 2D or 3D registration
#' is applied, and source images with one more dimension than the target (i.e.
#' 4D to 3D, or 3D to 2D) will be registered volumewise or slicewise, as
#' appropriate. In the latter case the last dimension of the resulting image is
#' taken from the source image, while all other dimensions come from the
#' target. One image of control points is returned for each registration
#' performed.
#'
#' @inheritParams niftyreg 
#' @param nLevels A single integer specifying the number of levels of the
#'   algorithm that should be applied. If zero, no optimisation will be
#'   performed, and the final control-point image will be the same as its
#'   initialisation value.
#' @param maxIterations A single integer specifying the maximum number of
#'   iterations to be used within each level. Fewer iterations may be used if a
#'   convergence test deems the process to have completed.
#' @param nBins A single integer giving the number of bins to use for the joint
#'   histogram created by the algorithm.
#' @param bendingEnergyWeight A numeric value giving the weight of the bending
#'   energy term in the cost function.
#' @param linearEnergyWeight A numeric value giving the weight of the linear
#'   energy term in the cost function.
#' @param jacobianWeight A numeric value giving the weight of the Jacobian
#'   determinant term in the cost function.
#' @param finalSpacing A numeric vector giving the spacing of control points in
#'   the final grid, along the X, Y and Z directions respectively. This is set
#'   from the initial control point image, if one is supplied.
#' @param spacingUnit A character string giving the units in which the
#'   \code{finalSpacing} is specified: either \code{"voxel"} for pixels/voxels,
#'   or \code{"world"} for real-world units (see \code{\link{pixunits}}).
#' @param verbose A single logical value: if \code{TRUE}, the code will give
#'   some feedback on its progress; otherwise, nothing will be output while the
#'   algorithm runs. Run time can be seconds or more, depending on the size and
#'   dimensionality of the images.
#' @return See \code{\link{niftyreg}}.
#' 
#' @note Performing a linear registration first, and then initialising the
#'   nonlinear transformation with the result (via the \code{init} parameter),
#'   is highly recommended in most circumstances.
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg}}, which can be used as an interface to this
#'   function, and \code{\link{niftyreg.linear}} for linear registration. Also,
#'   \code{\link{forward}} and \code{\link{reverse}} to extract
#'   transformations, and \code{\link{applyTransform}} to apply them to new
#'   images or points.
#' @references The algorithm used by this function is described in the
#' following publication.
#' 
#' M. Modat, G.R. Ridgway, Z.A. Taylor, M. Lehmann, J. Barnes, D.J. Hawkes,
#' N.C. Fox & S. Ourselin (2010). Fast free-form deformation using graphics
#' processing units. Computer Methods and Programs in Biomedicine
#' 98(3):278-284.
#' @export
niftyreg.nonlinear <- function (source, target, init = NULL, sourceMask = NULL, targetMask = NULL, symmetric = TRUE, nLevels = 3L, maxIterations = 150L, nBins = 64L, bendingEnergyWeight = 0.001, linearEnergyWeight = 0.01, jacobianWeight = 0, finalSpacing = c(5,5,5), spacingUnit = c("voxel","world"), interpolation = 3L, verbose = FALSE, estimateOnly = FALSE, sequentialInit = FALSE, internal = FALSE)
{
    if (missing(source) || missing(target))
        stop("Source and target images must be given")
    
    source <- .Call("retrieveImage", source, PACKAGE="RNiftyReg")
    target <- .Call("retrieveImage", target, PACKAGE="RNiftyReg")
    nSourceDim <- ndim(source)
    nTargetDim <- ndim(target)
    
    if (any(c(bendingEnergyWeight,jacobianWeight) < 0))
        stop("Penalty term weights must be nonnegative")
    if (bendingEnergyWeight + jacobianWeight > 1)
        stop("Penalty term weights cannot add up to more than 1")
    if (!(interpolation %in% c(0,1,3)))
        stop("Final interpolation specifier must be 0, 1 or 3")
    
    if (nLevels == 0)
        symmetric <- FALSE
    
    nReps <- ifelse(nSourceDim > nTargetDim, dim(source)[nSourceDim], 1L)
    spacingUnit <- match.arg(spacingUnit)
    spacingChanged <- FALSE
    
    if (!is.list(init))
        init <- list(init)
    if (length(init) != nReps)
    {
        if (sequentialInit)
            init <- c(init, rep(list(NULL),nReps-length(init)))
        else
            init <- rep(init, length.out=nReps)
    }
    init <- lapply(init, function(x) {
        if (is.null(x))
            return (x)
        else if (isImage(x))
        {
            currentSpacing <- pixdim(x)[1:3] / 2^max(0,nLevels-1)
            if (spacingChanged && !isTRUE(all.equal(currentSpacing, finalSpacing)))
                stop("Initial control point images must all use the same grid")
            finalSpacing <<- currentSpacing
            spacingUnit <<- "mm"
            spacingChanged <<- TRUE
            return (x)
        }
        else if (!isAffine(x))
            stop("Initial transform should be a control point image or affine matrix")
        else
            return (x)
    })
    
    if (spacingUnit == "voxel")
    {
        indices <- 1:min(3,nTargetDim)
        finalSpacing[indices] <- finalSpacing[indices] * abs(pixdim(target)[indices])
    }
    
    if (nTargetDim == 2)
        finalSpacing <- c(finalSpacing[1:2], 1)
    else
        finalSpacing <- finalSpacing[1:3]
    
    result <- .Call("regNonlinear", source, target, symmetric, nLevels, maxIterations, interpolation, sourceMask, targetMask, init, nBins, finalSpacing, bendingEnergyWeight, linearEnergyWeight, jacobianWeight, verbose, estimateOnly, sequentialInit, internal, PACKAGE="RNiftyReg")
    class(result) <- "niftyreg"
    
    return (result)
}


#' @rdname niftyreg
#' @export
as.array.niftyreg <- function (x, ...)
{
    as.array(x$image)
}


#' Extract forward and reverse transformations
#' 
#' These functions extract forward and reverse transformations in a form
#' compatible with \code{\link{applyTransform}} and other functions. They are
#' (S3) generic, but only methods for \code{"niftyreg"} objects currently
#' exist.
#' 
#' @param object An R object.
#' @param i The transformation number to extract. There will only be more than
#'   one in the case of multiple registration.
#' @param ... Additional arguments. Not currently used.
#' @return A transformation object, an image or affine matrix, with suitable
#'   attributes giving pointers to source and target images. If there is no
#'   transformation information in the object then \code{NULL} is returned.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg}}, \code{\link{applyTransform}}
#' @export
forward <- function (object, ...)
{
    UseMethod("forward")
}


#' @rdname forward
#' @export
forward.niftyreg <- function (object, i = 1, ...)
{
    if (is.null(object$forwardTransforms))
        return (NULL)
    else
    {
        result <- object$forwardTransforms[[i]]
        attr(result, "source") <- object$source[[i]]
        attr(result, "target") <- object$target
        return (result)
    }
}


#' @rdname forward
#' @export
reverse <- function (object, ...)
{
    UseMethod("reverse")
}


#' @rdname forward
#' @export
reverse.niftyreg <- function (object, i = 1, ...)
{
    if (is.null(object$reverseTransforms))
        return (NULL)
    else
    {
        result <- object$reverseTransforms[[i]]
        attr(result, "source") <- object$target
        attr(result, "target") <- object$source[[i]]
        return (result)
    }
}


#' Similarity measures between images
#' 
#' This function calculates a similarity measure between two images, after
#' resampling one into the space of the other. The only supported measure is
#' currently normalised mutual information, which is also used as a cost
#' function by the registration algorithms.
#' 
#' @param source The source image, in any acceptable form.
#' @param target The target image. Must have the same dimensionality as the
#'   source image.
#' @param targetMask An optional mask image in target space, whose nonzero
#'   region will be the area over which the measure is calculated.
#' @param interpolation A single integer specifying the type of interpolation
#'   to be applied to the source image when resampling it into the space of the
#'   target image. May be 0 (nearest neighbour), 1 (trilinear) or 3 (cubic
#'   spline). No other values are valid.
#' @return A single numeric value representing the similarity between the
#'   images.
#' 
#' @author Jon Clayden <code@@clayden.org>
#' @seealso \code{\link{niftyreg}}
#' @export
similarity <- function (source, target, targetMask = NULL, interpolation = 3L)
{
    if (!(interpolation %in% c(0,1,3)))
        stop("Final interpolation specifier must be 0, 1 or 3")
    
    return (.Call("calculateMeasure", source, target, targetMask, interpolation, PACKAGE="RNiftyReg"))
}
