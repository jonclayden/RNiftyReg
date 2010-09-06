\name{niftyreg}
\alias{niftyreg}

\title{Three dimensional image registration}
\description{
The \code{niftyreg} function performs registration for three dimensional images. Rigid-body (6 degrees of freedom) and affine (12 degrees of freedom) registration can currently be performed.
}

\usage{
niftyreg(source, target, targetMask = NULL, initAffine = NULL, scope = c("affine", "rigid"), nLevels = 3, maxIterations = 5, useBlockPercentage = 50, verbose = FALSE)
}

\arguments{
  \item{source}{The source image, an object of class \code{"nifti"}. Package \code{oro.nifti} defines this class and provides functions for reading and writing NIfTI files.}
  \item{target}{The target image, also an object of class \code{"nifti"}.}
  \item{targetMask}{An optional mask image (again a \code{"nifti"} object), whose nonzero region will be taken as the region of interest for the registration. Must have the same voxel and image dimensions as the target image.}
  \item{initAffine}{An optional affine matrix to initialise the algorithm. If \code{NULL}, the identity matrix is used, with an appropriate offset to account for differences in the image origins.}
  \item{scope}{A string describing the scope, or number of degrees of freedom (DOF), of the registration. Only \code{"affine"} (12 DOF) and \code{"rigid"} (6 DOF) are currently supported.}
  \item{nLevels}{A single integer specifying the number of levels of the algorithm that should be applied. If zero, no optimisation will be performed, and the final affine matrix will be the same as its initialisation value.}
  \item{maxIterations}{A single integer specifying the maximum number of iterations to be used within each level. Fewer iterations may be used if a convergence test deems the process to have completed.}
  \item{useBlockPercentage}{A single integer giving the percentage of the block to use.}
  \item{verbose}{A single logical value: if \code{TRUE}, the code will give some feedback on its progress; otherwise, nothing will be output while the algorithm runs.}
}

\details{
This function performs the dual operations of finding a transformation to optimise image alignment, and resampling the source image into the space of the target image.

The algorithm is based on a block-matching approach and Least Trimmed Squares (LTS) fitting. Firstly, the block matching provides a set of corresponding points between a target and a source image. Secondly, using this set of corresponding points, the best rigid or affine transformation is evaluated. This two-step loop is repeated until convergence to the best transformation.

In the NiftyReg implementation, normalised cross-correlation between the target and source blocks is used to evaluate correspondence. The block width is constant and has been set to 4 voxels. A coarse-to-fine approach is used, where the registration is first performed on down-sampled images (using a Gaussian filter to resample images), and finally performed on full resolution images.
}

\value{
A list of class \code{"niftyreg"} with components
\item{image}{An image object of class \code{"nifti"} representing the registered and resampled \code{source} image in the space of the \code{target} image.}
\item{affine}{A 4x4 matrix containing the optimised affine transformation matrix.}
\item{scope}{Copied from the function argument of the same name.}
}

\note{
A precalculated affine matrix can be applied to transform a new image by giving an appropriate \code{initAffine} and setting \code{nLevels} to 0.

If substantial parts of the target image are zero-valued, for example because the target image has been brain-extracted, it can be useful to pass it as a target mask as well as the target image. That is, \code{niftyreg(source, target, target)}.
}

\references{
The algorithm used by this function is described in the following publications.

Ourselin, S., Roche, A., Subsol, G., Pennec, X. and Ayache, N. (2000). Reconstructing a 3D structure from serial histological sections. Image and Vision Computing 19(1-2):25-31.

Ourselin, S., Stefanescu, R. and Pennec, X. (2002). Robust registration of multi-modal images: towards real-time clinical applications. Medical Image Computing and Computer-Assisted Intervention. Vol. 2489 of Lecture Notes in Computer Science, pp. 140-147.
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{nifti}} (no relation!), in the \code{oro.nifti} package, for creating the image objects passed to this function}
