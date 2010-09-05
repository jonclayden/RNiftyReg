\name{decomposeAffine}
\alias{decomposeAffine}

\title{Decompose an affine matrix into its constituent transformations}
\description{
An affine matrix is composed of translation, scale, skew and rotation transformations. This function extracts these components.
}

\usage{
decomposeAffine(affine, source = NULL, target = NULL, type = NULL)
}

\arguments{
  \item{affine}{A 4x4 matrix representing an affine transformation matrix.}
  \item{source}{A \code{"nifti"} object representing the source image for the transformation. Only required if \code{type} (or the matrix's \code{affineType} attribute) is \code{"niftyreg"}, in which case it is passed to \code{\link{convertAffine}}.}
  \item{target}{A \code{"nifti"} object representing the target image for the transformation. Only required if \code{type} (or the matrix's \code{affineType} attribute) is \code{"niftyreg"}, in which case it is passed to \code{\link{convertAffine}}.}
  \item{type}{The storage convention type of the affine matrix, if it is not stored in the \code{affineType} attribute of the matrix.}
}

\value{
A list with components
\item{scaleMatrix}{A 3x3 matrix representing only the scale operation embodied in the full affine transformation.}
\item{skewMatrix}{A 3x3 matrix representing only the skew operation embodied in the full affine transformation.}
\item{rotationMatrix}{A 3x3 matrix representing only the rotation operation embodied in the full affine transformation.}
\item{translation}{A length-3 named numeric vector representing the translations (generally in mm) in each of the X, Y and Z directions.}
\item{scales}{A length-3 named numeric vector representing the scale factors in each of the X, Y and Z directions. Scale factors of 1 represent no effect.}
\item{skews}{A length-3 named numeric vector representing the skews in each of the XY, XZ and YZ planes.}
\item{angles}{A length-3 named numeric vector representing the rotation angles (in radians) about each of the X, Y and Z directions; i.e. roll, pitch and yaw.}
}

\note{
The decomposition is not perfect, and there is one particular degenerate case when the pitch angle is very close to \code{pi/2} radians, known as ``Gimbal lock''. In this case the yaw angle is arbitrarily set to zero.

Affine matrices embodying rigid-body transformations include only 6 degrees of freedom, rather than the full 12, so skews will always be zero and scales will always be unity (to within rounding error).
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{convertAffine}}, \code{\link{niftyreg}}}
