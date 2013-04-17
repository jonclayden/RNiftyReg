\name{convertAffine}
\alias{convertAffine}

\title{Convert an affine matrix between different storage conventions}
\description{
Affine transformation matrices can be stored using a number of different conventions. This function converts between conventions.
}

\usage{
convertAffine(affine, source, target, newType = c("niftyreg","fsl"),
    currentType = NULL)
}

\arguments{
  \item{affine}{A 4x4 matrix representing an affine transformation matrix.}
  \item{source}{A \code{"nifti"} object representing the source image for the transformation.}
  \item{target}{A \code{"nifti"} object representing the target image for the transformation.}
  \item{newType}{The convention type to convert to. Currently valid values are \code{"niftyreg"} and \code{"fsl"} (for FSL FLIRT).}
  \item{currentType}{The current type of the affine matrix, which can be either of the valid values for \code{newType}. If \code{NULL}, the code attempts to determine the current type from the \code{affineType} attribute of the matrix.}
}

\details{
The source and target images for the original transformation are required to find the corresponding ``xform'' matrices, which are needed to perform the conversion.
}

\value{A copy of the \code{affine} object provided, converted to the requested storage convention.}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{readAffine}}, \code{\link{writeAffine}}}
