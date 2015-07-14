xform <- function (image, useQuaternionFirst = TRUE)
{
    return (.Call("getXform", image, isTRUE(useQuaternionFirst), PACKAGE="RNiftyReg"))
}

voxelToWorld <- function (points, image, simple = FALSE, ...)
{
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        voxelDims <- pixdim(image)[seq_len(ncol(points))]
        return (drop(t(apply(points-1, 1, function(x) x*abs(voxelDims)))))
    }
    else
    {
        affine <- xform(image, ...)
        return (applyAffine(points-1, affine))
    }
}

worldToVoxel <- function (points, image, simple = FALSE, ...)
{
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        voxelDims <- pixdim(image)[seq_len(ncol(points))]
        return (drop(t(apply(points, 1, function(x) x/abs(voxelDims)) + 1)))
    }
    else
    {
        affine <- solve(xform(image, ...))
        return (applyAffine(points, affine) + 1)
    }
}
