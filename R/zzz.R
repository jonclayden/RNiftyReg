#' @import ore RNifti
#' @importFrom Rcpp evalCpp
#' @importFrom splines bs
#' @importFrom stats lm na.omit predict
#' @importFrom utils read.table
#' @useDynLib RNiftyReg, .registration = TRUE, .fixes = "C_"
NULL

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
