\name{writeAffine}
\alias{writeAffine}

\title{Write an affine matrix to a file}
\description{
This function is used to write a 4x4 numeric matrix representing an affine transformation to a file. The \code{affineType} attribute is also written, as a comment.
}

\usage{
writeAffine(affine, fileName)
}

\arguments{
  \item{affine}{A 4x4 affine matrix.}
  \item{fileName}{A string giving the file name to write the matrix to.}
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{write.table}}, \code{\link{readAffine}}, \code{\link{convertAffine}}}
