\name{transformWithControlPoints}
\alias{transformWithControlPoints}

\title{Transform points using a control point image}
\description{
This function is used to transform points from source to target space using a control point image, possibly obtained from \code{\link{niftyreg.nonlinear}}. A subvoxel location is target space is returned. 
}

\usage{
transformWithControlPoints(points, controlPointImage, source, target,
    nearest = FALSE)
}

\arguments{
  \item{points}{A vector giving the coordinates of a point, or matrix with one point per row. These must be given as voxel-based locations: see \code{\link{transformWorldToVoxel}} for a function to convert from ``real world'' coordinates.}
  \item{controlPointImage}{An object of class \code{"nifti"} containing the control points parameterising the transformation.}
  \item{source}{A \code{"nifti"} object representing the source image for the transformation.}
  \item{target}{A \code{"nifti"} object representing the target image for the transformation.}
  \item{nearest}{A logical value: if \code{TRUE}, the return value contains just the location of the nearest voxel to each requested location. This is computationally more efficient than simply rounding the result.}
}

\details{
The transformation performed by this function is not exact, and its accuracy will depend to some extent on the density of the control point grid and the geometry of the deformation in the vicinity of the points of interest. Nevertheless, it should be quite sufficient for most purposes.

The method is to first convert the control points to a deformation field, which encodes the location of each target space voxel in the source space. The target voxel closest to the requested location is found by searching through this deformation field, and returned if \code{nearest} is \code{TRUE} or it coincides exactly with the requested location. Otherwise, a block of four voxels in each dimension around the point of interest is extracted from the deformation field, and the final location is estimated by cubic spline regression.
}

\value{A vector or matrix of transformed points.}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{niftyreg.nonlinear}}}
