\name{transformWithAffine}
\alias{transformWithAffine}

\title{Transform points using an affine matrix}
\description{
This function is used to transform points from source to target space using an affine transformation, possibly obtained from \code{\link{niftyreg.linear}}. The exact subvoxel location in target space is returned.
}

\usage{
transformWithAffine(points, affine, source, target, type = NULL)
}

\arguments{
  \item{points}{A vector giving the coordinates of a point, or a matrix with one point per row. These must be given as voxel-based locations: see \code{\link{transformWorldToVoxel}} for a function to convert from ``real world'' coordinates.}
  \item{affine}{A 4x4 matrix representing an affine transformation matrix.}
  \item{source}{A \code{"nifti"} object representing the source image for the transformation.}
  \item{target}{A \code{"nifti"} object representing the target image for the transformation.}
  \item{type}{The convention type of the affine matrix. Currently valid values are \code{"niftyreg"} and \code{"fsl"} (for FSL FLIRT). If \code{NULL}, the code attempts to determine the current type from the \code{affineType} attribute of the matrix.}
}

\value{A vector or matrix of transformed points.}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{For transforming points using the xform matrix stored in an image, \code{\link{transformVoxelToWorld}} provides a simpler interface.}
