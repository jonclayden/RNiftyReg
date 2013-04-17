\name{readAffine}
\alias{readAffine}

\title{Read an affine matrix from a file}
\description{
This function is used to read a 4x4 numeric matrix representing an affine transformation from a file. It is a wrapper around \code{read.table} which additionally ensures that the \code{affineType} attribute is set. This is required because there are different conventions for storing affine matrices.
}

\usage{
readAffine(fileName, type = NULL)
}

\arguments{
  \item{fileName}{A string giving the file name to read the affine matrix from.}
  \item{type}{The type of the affine matrix, which describes what convention is it is stored with. Currently valid values are \code{"niftyreg"} and \code{"fsl"} (for FSL FLIRT). If \code{NULL}, the function will look in the file for a comment specifying the type.}
}

\value{An affine matrix with \code{affineType} attribute set appropriately.}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{read.table}}, \code{\link{writeAffine}}, \code{\link{convertAffine}}}
