\name{getDeformationField}
\alias{getDeformationField}

\title{Calculate the deformation field for a transformation}
\description{
This function is used to calculate the deformation field corresponding to a specified linear or nonlinear transformation. The deformation field gives the location in source image space corresponding to the centre of each voxel in target space.
}

\usage{
getDeformationField(target, affine = NULL, controlPointImage = NULL,
    jacobian = TRUE)
}

\arguments{
  \item{target}{A \code{"nifti"} object representing the target image for the transformation.}
  \item{affine}{A 4x4 matrix representing an affine matrix, using the \code{"niftyreg"} convention. Either this or a control point image must be specified.}
  \item{controlPointImage}{A \code{"nifti"} object representing the control points for a nonlinear transformation. Either this or an affine matrix must be specified.}
  \item{jacobian}{A logical value: if \code{TRUE}, a Jacobian determinant map is also calculated and returned.}
}

\value{
A list with components
\item{deformationField}{An object of class \code{"nifti"} representing the deformation field.}
\item{jacobian}{An object of class \code{"nifti"} representing the Jacobian determinant map, if requested. Otherwise this element is not present.}
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{niftyreg.linear}}, \code{\link{niftyreg.nonlinear}}}
