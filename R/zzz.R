#' @import ore RNifti
#' @importFrom Rcpp evalCpp
#' @importFrom splines bs
#' @importFrom stats lm na.omit predict
#' @importFrom utils object.size read.table
#' @useDynLib RNiftyReg
.onLoad <- function (libname, pkgname)
{
    .Call("initialise", PACKAGE="RNiftyReg")
}
