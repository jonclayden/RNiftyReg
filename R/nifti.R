readNifti <- function (file, internal = FALSE)
{
    if (!is.character(file))
        report(OL$Error, "File name(s) must be specified in a character vector")
    if (length(file) == 0)
        report(OL$Error, "File name vector is empty")
    else if (length(file) > 1)
        lapply(file, function(f) .Call("readNifti", f, internal, PACKAGE="RNiftyReg"))
    else
        .Call("readNifti", file, internal, PACKAGE="RNiftyReg")
}

writeNifti <- function (image, file)
{
    .Call("writeNifti", image, file, PACKAGE="RNiftyReg")
}
