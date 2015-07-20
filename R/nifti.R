readNifti <- function (file, internal = FALSE)
{
    if (!is.character(file))
        stop("File name(s) must be specified in a character vector")
    if (length(file) == 0)
        stop("File name vector is empty")
    else if (length(file) > 1)
        lapply(file, function(f) .Call("readNifti", f, internal, PACKAGE="RNiftyReg"))
    else
        .Call("readNifti", file, internal, PACKAGE="RNiftyReg")
}

writeNifti <- function (image, file)
{
    .Call("writeNifti", image, file, PACKAGE="RNiftyReg")
}
