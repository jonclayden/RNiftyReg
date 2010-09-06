\name{xformToAffine}
\alias{xformToAffine}

\title{Create an affine matrix corresponding to the ``xform'' of an image}
\description{
This function converts the ``qform'' or ``sform'' information in a NIfTI header into its corresponding affine matrix. These two ``xform'' mechanisms are defined by the NIfTI standard and may both be in use in a particular image header.
}

\usage{
xformToAffine(image, useQuaternionFirst = TRUE)
}

\arguments{
  \item{image}{A object of class \code{"nifti"} from which to take the information.}
  \item{useQuaternionFirst}{A single logical value. If \code{TRUE}, the ``qform'' matrix will be used first, if it is defined; otherwise the ``sform'' matrix will take priority.}
}

\value{A affine matrix corresponding to the ``qform'' or ``sform'' information in the image header.}

\references{The NIfTI-1 standard (\url{http://nifti.nimh.nih.gov/nifti-1}) is the definitive reference on ``xform'' conventions.}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{nifti}, in the \code{oro.nifti} package}
