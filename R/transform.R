transformWithAffine <- function (points, affine, source = NULL, target = NULL, type = NULL)
{
    affine <- convertAffine(affine, source, target, "fsl", type)
    
    points <- transformVoxelToWorld(points, source, simple=TRUE)
    
    if (!is.matrix(points))
        points <- matrix(points, nrow=1)
    
    nDims <- ncol(points)
    if (nDims != 2 && nDims != 3)
        report(OL$Error, "Points must be two or three dimensional")
    
    for (i in (nDims+1):4)
        points <- cbind(points, 1)
    newPoints <- affine %*% t(points)
    newPoints <- drop(t(newPoints[1:3,,drop=FALSE]))
    newPoints <- transformWorldToVoxel(newPoints, target, simple=TRUE)
    
    return (newPoints)
}

transformVoxelToWorld <- function (points, image, simple = FALSE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    if (image@dim_[1] != 3)
        report(OL$Error, "Only three-dimensional images can be used at present")
    
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        return (drop(apply(points-1, 1, function(x) x*abs(image@pixdim[2:4]))))
    }
    else
    {
        affine <- xformToAffine(image, ...)
        return (transformWithAffine(points-1, affine, type="fsl"))
    }
}

transformWorldToVoxel <- function (points, image, simple = FALSE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    if (image@dim_[1] != 3)
        report(OL$Error, "Only three-dimensional images can be used at present")
    
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        return (drop(apply(points, 1, function(x) x/abs(image@pixdim[2:4])) + 1))
    }
    else
    {
        affine <- solve(xformToAffine(image, ...))
        return (transformWithAffine(points, affine, type="fsl") + 1)
    }
}

transformWithControlPoints <- function (points, controlPointImage, source = NULL, target = NULL, nearest = FALSE)
{
    if (!is.nifti(controlPointImage))
        report(OL$Error, "Control point image must be specified as a \"nifti\" object")
    
    points <- transformVoxelToWorld(points, source)
    
    if (!is.matrix(points))
        points <- matrix(points, nrow=1)
    
    nDims <- ncol(points)
    if (nDims != 2 && nDims != 3)
        report(OL$Error, "Points must be two or three dimensional")
    
    # This function takes world coordinates but returns voxels
    newPoints <- drop(.Call("cp_transform_R", .fixTypes(controlPointImage), .fixTypes(target), points, as.logical(nearest), PACKAGE="RNiftyReg"))
    
    return (newPoints)
}
