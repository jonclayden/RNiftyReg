.First.lib <- function (libname, pkgname)
{
    library.dynam("RNiftyReg", package="RNiftyReg")
    
    if (is.null(getOption("niftiAuditTrail")))
        options(niftiAuditTrail=FALSE)
}
