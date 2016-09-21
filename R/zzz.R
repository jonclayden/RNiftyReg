#' @import ore RNifti
#' @importFrom Rcpp evalCpp
#' @importFrom splines bs
#' @importFrom stats lm na.omit predict
#' @importFrom utils read.table
#' @export "pixdim<-" "pixunits<-" "qform<-" "sform<-" dumpNifti ndim pixdim pixunits retrieveNifti updateNifti voxelToWorld worldToVoxel writeNifti xform
#' @useDynLib RNiftyReg
.onLoad <- function (libname, pkgname)
{
    .Call("initialise", PACKAGE="RNiftyReg")
}
