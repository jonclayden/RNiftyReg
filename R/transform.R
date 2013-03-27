transformWithAffine <- function (points, affine, type)

transformVoxelToWorld <- function (points, image, useOrigin = TRUE, ...)
{
    affine <- xformToAffine(image, ...)
    if (!useOrigin)
        affine[1:3,4] <- 0
    
    return (solve(affine) %*% points)
}

transformWorldToVoxel <- function (points, image, useOrigin = TRUE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    if (image@dim_[1] != 3)
        report(OL$Error, "Only three-dimensional images can be used at present")
    if (!is.numeric(point) || length(point) != 3)
        report(OL$Error, "Point must be specified as a numeric vector of length 3")
    
    affine <- xformToAffine(image, ...)
    if (!useOrigin)
        affine[1:3,4] <- 0
    
    return (affine %*% points)
}
