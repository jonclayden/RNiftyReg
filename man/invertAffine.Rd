\name{invertAffine}
\alias{invertAffine}

\title{Invert an affine matrix}
\description{
This function is used to invert an affine matrix. It is a wrapper around \code{\link{solve}}, which additionally preserves the \code{affineType} attribute.
}

\usage{
invertAffine(affine)
}

\arguments{
  \item{affine}{An existing 4x4 affine matrix.}
}

\value{The inverted affine matrix.}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{solve}}}
