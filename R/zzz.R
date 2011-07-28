.onLoad <- function (libname, pkgname)
{
    if (is.null(getOption("niftiAuditTrail")))
        options(niftiAuditTrail=FALSE)
}
