\name{transformVoxelToWorld}
\alias{transformVoxelToWorld}
\alias{transformWorldToVoxel}

\title{Transform points between voxel and ``world'' coordinates}
\description{
These functions are used to transform points from dimensionless pixel or voxel coordinates to ``real-world'' coordinates, typically in millimetres, and back.
}

\usage{
transformVoxelToWorld(points, image, simple = FALSE, ...)

transformWorldToVoxel(points, image, simple = FALSE, ...)
}

\arguments{
  \item{points}{A vector giving the coordinates of a point, or matrix with one point per row.}
  \item{image}{The image in whose space the points are given, an object of class \code{"nifti"}.}
  \item{simple}{A logical value: if \code{TRUE} then the transformation is performed simply by rescaling the points according to the voxel dimensions recorded in the \code{image}. Otherwise the full xform matrix is used.}
  \item{...}{Additional arguments to \code{\link{xformToAffine}}.}
}

\value{A vector or matrix of transformed points.}

\note{
Voxel coordinates are assumed by these functions to use R's indexing convention, beginning from 1.
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{xformToAffine}}}
