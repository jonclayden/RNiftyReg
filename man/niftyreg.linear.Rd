\name{niftyreg.linear}
\alias{niftyreg.linear}
\alias{applyAffine}

\title{Two and three dimensional linear image registration}
\description{
The \code{niftyreg.linear} function performs linear registration for two and three dimensional images. 4D images may also be registered volumewise to a 3D image, or 3D images slicewise to a 2D image. Rigid-body (6 degrees of freedom) and affine (12 degrees of freedom) registration can currently be performed. A precalculated transformation can be applied to a new image using the \code{applyAffine} function.
}

\usage{
niftyreg.linear(source, target, targetMask = NULL, initAffine = NULL,
    scope = c("affine","rigid"), nLevels = 3, maxIterations = 5,
    useBlockPercentage = 50, finalInterpolation = 3, verbose = FALSE)

applyAffine(affine, source, target, affineType = NULL,
    finalInterpolation = 3)
}

\arguments{
  \item{source}{The source image, an object of class \code{"nifti"} with 2, 3 or 4 dimensions. Package \code{oro.nifti} defines this class and provides functions for reading and writing NIfTI files.}
  \item{target}{The target image, an object of class \code{"nifti"} with 2 or 3 dimensions.}
  \item{targetMask}{An optional mask image (again a \code{"nifti"} object), whose nonzero region will be taken as the region of interest for the registration. Must have the same voxel and image dimensions as the target image.}
  \item{initAffine}{An optional affine matrix, or list of matrices, to initialise the algorithm. If \code{NULL}, the identity matrix is used, with an appropriate offset to account for differences in the image origins.}
  \item{scope}{A string describing the scope, or number of degrees of freedom (DOF), of the registration. Only \code{"affine"} (12 DOF) and \code{"rigid"} (6 DOF) are currently supported.}
  \item{nLevels}{A single integer specifying the number of levels of the algorithm that should be applied. If zero, no optimisation will be performed, and the final affine matrix will be the same as its initialisation value.}
  \item{maxIterations}{A single integer specifying the maximum number of iterations to be used within each level. Fewer iterations may be used if a convergence test deems the process to have completed.}
  \item{useBlockPercentage}{A single integer giving the percentage of blocks to use for calculating correspondence at each step of the algorithm. The blocks with the highest intensity variance will be chosen.}
  \item{finalInterpolation}{A single integer specifying the type of interpolation to be applied to the final resampled image. May be 0 (nearest neighbour), 1 (trilinear) or 3 (cubic spline). No other values are valid.}
  \item{verbose}{A single logical value: if \code{TRUE}, the code will give some feedback on its progress; otherwise, nothing will be output while the algorithm runs.}
  \item{affine}{For \code{applyAffine}, the affine transformation(s) to apply to the source image.}
  \item{affineType}{For \code{applyAffine}, the storage convention type of the affine matrix, if it is not stored in the \code{affineType} attribute of the matrix.}
}

\details{
This function performs the dual operations of finding a transformation to optimise image alignment, and resampling the source image into the space of the target image.

The algorithm is based on a block-matching approach and Least Trimmed Squares (LTS) fitting. Firstly, the block matching provides a set of corresponding points between a target and a source image. Secondly, using this set of corresponding points, the best rigid or affine transformation is evaluated. This two-step loop is repeated until convergence to the best transformation.

In the NiftyReg implementation, normalised cross-correlation between the target and source blocks is used to evaluate correspondence. The block width is constant and has been set to 4 voxels. A coarse-to-fine approach is used, where the registration is first performed on down-sampled images (using a Gaussian filter to resample images), and finally performed on full resolution images.

The source image may have 2, 3 or 4 dimensions, and the target 2 or 3. The dimensionality of the target image determines whether 2D or 3D registration is applied, and source images with one more dimension than the target (i.e. 4D to 3D, or 3D to 2D) will be registered volumewise or slicewise, as appropriate. In the latter case the last dimension of the resulting image is taken from the source image, while all other dimensions come from the target. One affine matrix is returned for each registration performed.

The \code{applyAffine} function is a convenience wrapper that calls \code{niftyreg.linear} with \code{nLevels=0} to apply the specified transformation without any further optimisation. Note that a target image must still be specified in this case, since the metadata associated with that image is needed by \code{niftyreg.linear}.
}

\value{
See \code{\link{niftyreg}}.
}

\note{
If substantial parts of the target image are zero-valued, for example because the target image has been brain-extracted, it can be useful to pass it as a target mask as well as the target image, viz. \code{niftyreg.linear(source, target, target)}.

There is no reason that arrays that do not represent medical images cannot be registered using this function. A standard R array can be converted to a valid \code{"nifti"} object easily for these purposes using the \code{as.nifti} function in the \code{oro.nifti} package.
}

\references{
The algorithm used by this function is described in the following publications.

S. Ourselin, A. Roche, G. Subsol, X. Pennec & N. Ayache (2000). Reconstructing a 3D structure from serial histological sections. Image and Vision Computing 19(1-2):25-31.

S. Ourselin, R. Stefanescu & X. Pennec (2002). Robust registration of multi-modal images: towards real-time clinical applications. Medical Image Computing and Computer-Assisted Intervention. Vol. 2489 of Lecture Notes in Computer Science, pp. 140-147.
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{niftyreg}}, which can be used as an interface to this function, and \code{\link{niftyreg.nonlinear}} for nonlinear registration. See \code{nifti} (no relation!), in the \code{oro.nifti} package, for creating the image objects passed to this function. Useful related functions are \code{as.nifti}, \code{readNIfTI} and \code{writeNIfTI}.}
