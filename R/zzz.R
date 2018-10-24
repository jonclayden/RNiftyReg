#' @import ore RNifti
#' @importFrom Rcpp evalCpp
#' @importFrom splines bs
#' @importFrom stats lm na.omit predict
#' @importFrom utils read.table packageVersion
#' @useDynLib RNiftyReg, .registration = TRUE, .fixes = "C_"
.onLoad <- function (libname, pkgname)
{
    # RNifti v0.10.0 changed the fields in a NiftiImage C++ object, which means old external pointers are invalid and likely to cause crashes
    buildVersion <- .Call(C_RNifti_version)
    currentVersion <- utils::packageVersion("RNifti")
    if ((buildVersion < 10) != (currentVersion < "0.10.0"))
        warning("This package was built against an RNifti API version that does not match the current one - a clean reinstall of RNiftyReg is advisable", call.=FALSE)
}

#' @export
RNifti::`pixdim<-`
#' @export
RNifti::`pixunits<-`
#' @export
RNifti::`qform<-`
#' @export
RNifti::`sform<-`
#' @export
RNifti::dumpNifti
#' @export
RNifti::ndim
#' @export
RNifti::pixdim
#' @export
RNifti::pixunits
#' @export
RNifti::retrieveNifti
#' @export
RNifti::updateNifti
#' @export
RNifti::voxelToWorld
#' @export
RNifti::worldToVoxel
#' @export
RNifti::writeNifti
#' @export
RNifti::xform
