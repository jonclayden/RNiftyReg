\name{niftyreg.nonlinear}
\alias{niftyreg.nonlinear}
\alias{applyControlPoints}

\title{Two and three dimensional nonlinear image registration}
\description{
The \code{niftyreg.nonlinear} function performs nonlinear registration for two and three dimensional images. 4D images may also be registered volumewise to a 3D image, or 3D images slicewise to a 2D image. The warping is based on free-form deformations, parameterised using an image of control points. A precalculated transformation can be applied to a new image using the \code{applyControlPoints} function.
}

\usage{
niftyreg.nonlinear(source, target, targetMask = NULL, initAffine = NULL,
    initControl = NULL, symmetric = FALSE, sourceMask = NULL, nLevels = 3,
    maxIterations = 300, nBins = 64, bendingEnergyWeight = 0.005,
    jacobianWeight = 0, inverseConsistencyWeight = 0.01,
    finalSpacing = c(5,5,5), spacingUnit = c("vox","mm"),
    finalInterpolation = 3, verbose = FALSE, estimateOnly = FALSE)

applyControlPoints(controlPointImage, source, target,
    finalInterpolation = 3)
}

\arguments{
  \item{source}{The source image, an object of class \code{"nifti"} with 2, 3 or 4 dimensions. Package \code{oro.nifti} defines this class and provides functions for reading and writing NIfTI files.}
  \item{target}{The target image, an object of class \code{"nifti"} with 2 or 3 dimensions.}
  \item{targetMask}{An optional mask image (again a \code{"nifti"} object), whose nonzero region will be taken as the region of interest for the registration. Must have the same voxel and image dimensions as the target image.}
  \item{initAffine}{An optional affine matrix, or list of matrices, to initialise the algorithm. If both this parameter and \code{initControl} are \code{NULL}, the identity matrix is used, with an appropriate offset to account for differences in the image origins.}
  \item{initControl}{An optional image of class \code{"nifti"}, or a list of images, representing fields of previously-calculated control points to use as initialisation for the algorithm. This parameter takes priority over \code{initAffine} if both are not \code{NULL}, but may not be specified if \code{symmetric} is \code{TRUE}.}
  \item{symmetric}{A single logical value: if \code{TRUE}, a symmetric version of the registration algorithm is used to simultaneously obtain transformations from \code{source} to \code{target} and \code{target} to \code{source}. There are some restrictions in this mode: see Details.}
  \item{sourceMask}{An optional mask image whose nonzero region will be taken as the region of interest for the registration. Must have the same voxel and image dimensions as the source image. Ignored if \code{symmetric} is \code{FALSE}.}
  \item{nLevels}{A single integer specifying the number of levels of the algorithm that should be applied. If zero, no optimisation will be performed, and the final control-point image will be the same as its initialisation value.}
  \item{maxIterations}{A single integer specifying the maximum number of iterations to be used within each level. Fewer iterations may be used if a convergence test deems the process to have completed.}
  \item{nBins}{A single integer giving the number of bins to use for the joint histogram created by the algorithm.}
  \item{bendingEnergyWeight}{A numeric value giving the weight of the bending energy term in the cost function.}
  \item{jacobianWeight}{A numeric value giving the weight of the Jacobian determinant term in the cost function.}
  \item{inverseConsistencyWeight}{A numeric value giving the weight of the term ensuring inverse consistency in the cost function. Ignored if \code{symmetric} is \code{FALSE}.}
  \item{finalSpacing}{A numeric vector giving the spacing of control points in the final grid, along the X, Y and Z directions respectively. This is set from the \code{initControl} image if one is supplied.}
  \item{spacingUnit}{A character string giving the units in which the \code{finalSpacing} is specified: either \code{"vox"} for voxels or \code{"mm"} for millimetres (which is assumed to be the spatial unit of the source and target images).}
  \item{finalInterpolation}{A single integer specifying the type of interpolation to be applied to the final resampled image. May be 0 (nearest neighbour), 1 (trilinear) or 3 (cubic spline). No other values are valid.}
  \item{verbose}{A single logical value: if \code{TRUE}, the code will give some feedback on its progress; otherwise, nothing will be output while the algorithm runs.}
  \item{estimateOnly}{A single logical value: if \code{TRUE}, the transformation(s) will be estimated but the image(s) will not be resampled.}
  \item{controlPointImage}{For \code{applyControlPoints}, the control point map to apply to the source image.}
}

\details{
This function performs the dual operations of finding a transformation to optimise image alignment, and resampling the source image into the space of the target image (and vice-versa, if \code{symmetric} is \code{TRUE}). Unlike \code{\link{niftyreg.linear}}, this transformation is nonlinear, and the degree of deformation may vary across the image.

The nonlinear warping is based on free-form deformations. A lattice of equally-spaced control points is defined over the target image, each of which can be moved to locally modify the mapping to the source image. In order to assess the quality of the warping between the two images, an objective function based on the normalised mutual information is used, with penalty terms based on the bending energy or the squared log of the Jacobian determinant. The objective function value is optimised using a conjugate gradient scheme.

The source image may have 2, 3 or 4 dimensions, and the target 2 or 3. The dimensionality of the target image determines whether 2D or 3D registration is applied, and source images with one more dimension than the target (i.e. 4D to 3D, or 3D to 2D) will be registered volumewise or slicewise, as appropriate. In the latter case the last dimension of the resulting image is taken from the source image, while all other dimensions come from the target. One image of control points is returned for each registration performed.

The \code{symmetric} option allows registration to be performed in both directions simultaneously. This may be preferable to performing two one-way registrations, since a constraint is applied to ensure that the control point maps in the two directions are consistent. However, the source and target images must have the same dimensionality in this case, and an initial control point map may not be used.

The \code{applyControlPoints} function is a convenience wrapper that calls \code{niftyreg.nonlinear} with \code{nLevels=0} to apply the specified transformation without any further optimisation. Note that a target image must still be specified in this case, since the metadata associated with that image is needed by \code{niftyreg.nonlinear}.
}

\value{
See \code{\link{niftyreg}}.
}

\note{
Performing a linear registration first, and then initialising the nonlinear transformation with the result (via the \code{initAffine} parameter), is highly recommended in most circumstances.

If substantial parts of the target image are zero-valued, for example because the target image has been brain-extracted, it can be useful to pass it as a target mask as well as the target image, viz. \code{niftyreg.nonlinear(source, target, target)}.

There is no reason that arrays that do not represent medical images cannot be registered using this function. A standard R array can be converted to a valid \code{"nifti"} object easily for these purposes using the \code{as.nifti} function in the \code{oro.nifti} package.
}

\references{
The algorithm used by this function is described in the following publication.

M. Modat, G.R. Ridgway, Z.A. Taylor, M. Lehmann, J. Barnes, D.J. Hawkes, N.C. Fox & S. Ourselin (2010). Fast free-form deformation using graphics processing units. Computer Methods and Programs in Biomedicine 98(3):278-284.
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{niftyreg}}, which can be used as an interface to this function, and \code{\link{niftyreg.linear}} for linear registration. See \code{nifti} (no relation!), in the \code{oro.nifti} package, for creating the image objects passed to this function. Useful related functions are \code{as.nifti}, \code{readNIfTI} and \code{writeNIfTI}.}
