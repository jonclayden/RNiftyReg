\name{niftyreg}
\alias{niftyreg}

\title{Two and three dimensional image registration}
\description{
The \code{niftyreg} function performs linear or nonlinear registration for two and three dimensional images. 4D images may also be registered volumewise to a 3D image, or 3D images slicewise to a 2D image. This function is a common wrapper for \code{\link{niftyreg.linear}} and \code{\link{niftyreg.nonlinear}}.
}

\usage{
niftyreg(source, target, targetMask = NULL, initAffine = NULL,
            scope = c("affine","rigid","nonlinear"), ...)
}

\arguments{
  \item{source}{The source image, an object of class \code{"nifti"} with 2, 3 or 4 dimensions. Package \code{oro.nifti} defines this class and provides functions for reading and writing NIfTI files.}
  \item{target}{The target image, an object of class \code{"nifti"} with 2 or 3 dimensions.}
  \item{targetMask}{An optional mask image (again a \code{"nifti"} object), whose nonzero region will be taken as the region of interest for the registration. Must have the same voxel and image dimensions as the target image.}
  \item{initAffine}{An optional affine matrix to initialise the algorithm.}
  \item{scope}{A string describing the scope, or number of degrees of freedom (DOF), of the registration. The currently-supported values are \code{"affine"} (12 DOF), \code{"rigid"} (6 DOF) or \code{"nonlinear"} (high DOF, with the exact number depending on the image sizes).}
  \item{...}{Further arguments to \code{\link{niftyreg.linear}} or \code{\link{niftyreg.nonlinear}}.}
}

\value{
A list of class \code{"niftyreg"} with components
\item{image}{An image object of class \code{"nifti"} representing the registered and resampled \code{source} image in the space of the \code{target} image.}
\item{affine}{A list of 4x4 matrices containing the optimised affine transformations for each slice or volume of the source image. If the target is a 2D image, elements involving the Z dimension will have no effect. This element is \code{NULL} if nonlinear registration is performed.}
\item{control}{A list of objects of class \code{"nifti"}, representing the control point images for each warping. This element is \code{NULL} if linear (rigid or affine) registration is performed.}
\item{iterations}{A list of integer vectors giving the number of iterations actually run within each level, for each slice or volume of the source image. Note that for the first level of the linear algorithm specifically, twice the specified number of iterations is allowed.}
\item{scope}{Copied from the function argument of the same name.}
}

\note{
If substantial parts of the target image are zero-valued, for example because the target image has been brain-extracted, it can be useful to pass it as a target mask as well as the target image, viz. \code{niftyreg(source, target, target)}.

There is no reason that arrays that do not represent medical images cannot be registered using this function. A standard R array can be converted to a valid \code{"nifti"} object easily for these purposes using the \code{as.nifti} function in the \code{oro.nifti} package.
}

\references{
Please see \code{\link{niftyreg.linear}} or \code{\link{niftyreg.nonlinear}} for references relating to each type of registration.
}

\author{Jon Clayden <jon.clayden+rniftyreg@gmail.com>}

\seealso{\code{\link{niftyreg.linear}} and \code{\link{niftyreg.nonlinear}}, which do most of the work. See \code{nifti} (no relation!), in the \code{oro.nifti} package, for creating the image objects passed to this function. Useful related functions are \code{as.nifti}, \code{readNIfTI} and \code{writeNIfTI}.}
