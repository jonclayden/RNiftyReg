transformWithAffine <- function (points, affine, voxel = FALSE, source = NULL, target = NULL, type = NULL)
{
    if (!is.matrix(affine) || !isTRUE(all.equal(dim(affine), c(4,4))))
        report(OL$Error, "Specified affine matrix is not valid")
    if (is.null(type))
    {
        type <- attr(affine, "affineType")
        if (is.null(type))
            report(OL$Error, "The current affine type was not specified and is not stored with the matrix")
    }
    
    if (type == "niftyreg")
        affine <- convertAffine(affine, source, target, "fsl", "niftyreg")
    
    if (voxel)
        points <- transformVoxelToWorld(points, source, useOrigin=FALSE)
    
    if (!is.matrix(points))
        points <- matrix(points, nrow=1)
    
    nDims <- ncol(points)
    if (nDims != 2 && nDims != 3)
        report(OL$Error, "Points must be two or three dimensional")
    
    for (i in (nDims+1):4)
        points <- cbind(points, 1)
    newPoints <- affine %*% t(points)
    newPoints <- drop(t(newPoints[1:3,,drop=FALSE]))
    
    if (voxel)
        newPoints <- transformWorldToVoxel(newPoints, target, useOrigin=FALSE)
    
    return (newPoints)
}

transformVoxelToWorld <- function (points, image, useOrigin = TRUE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    if (image@dim_[1] != 3)
        report(OL$Error, "Only three-dimensional images can be used at present")
    
    affine <- xformToAffine(image, keepOrigin=useOrigin, ...)
    return (transformWithAffine(points-1, affine))
}

transformWorldToVoxel <- function (points, image, useOrigin = TRUE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    if (image@dim_[1] != 3)
        report(OL$Error, "Only three-dimensional images can be used at present")
    
    affine <- solve(xformToAffine(image, keepOrigin=useOrigin, ...))
    return (transformWithAffine(points, affine) + 1)
}
