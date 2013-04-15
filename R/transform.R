.applyAffine <- function (points, affine)
{
    if (!is.matrix(affine) || !isTRUE(all.equal(dim(affine), c(4,4))))
        report(OL$Error, "Specified affine matrix is not valid")
    
    if (!is.matrix(points))
        points <- matrix(points, nrow=1)
    
    nDims <- ncol(points)
    if (nDims != 2 && nDims != 3)
        report(OL$Error, "Points must be two or three dimensional")
    
    if (nDims == 2)
        affine <- matrix(affine[c(1,2,4,5,6,8,13,14,16)], ncol=3, nrow=3)
    
    points <- cbind(points, 1)
    newPoints <- affine %*% t(points)
    newPoints <- drop(t(newPoints[1:nDims,,drop=FALSE]))
    
    return (newPoints)
}

transformWithAffine <- function (points, affine, source = NULL, target = NULL, type = NULL)
{
    affine <- convertAffine(affine, source, target, "fsl", type)
    points <- transformVoxelToWorld(points, source, simple=TRUE)
    newPoints <- .applyAffine(points, affine)
    newPoints <- transformWorldToVoxel(newPoints, target, simple=TRUE)
    
    return (newPoints)
}

transformVoxelToWorld <- function (points, image, simple = FALSE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        voxelDims <- image@pixdim[seq_len(ncol(points))+1]
        return (drop(apply(points-1, 1, function(x) x*abs(voxelDims))))
    }
    else
    {
        affine <- xformToAffine(image, ...)
        return (.applyAffine(points-1, affine))
    }
}

transformWorldToVoxel <- function (points, image, simple = FALSE, ...)
{
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    
    if (simple)
    {
        if (!is.matrix(points))
            points <- matrix(points, nrow=1)
        voxelDims <- image@pixdim[seq_len(ncol(points))+1]
        return (drop(apply(points, 1, function(x) x/abs(voxelDims)) + 1))
    }
    else
    {
        affine <- solve(xformToAffine(image, ...))
        return (.applyAffine(points, affine) + 1)
    }
}

transformWithControlPoints <- function (points, controlPointImage, source = NULL, target = NULL, nearest = FALSE)
{
    library("splines")
    
    if (!is.nifti(controlPointImage))
        report(OL$Error, "Control point image must be specified as a \"nifti\" object")
    
    points <- transformVoxelToWorld(points, source)
    
    if (!is.matrix(points))
        points <- matrix(points, nrow=1)
    
    nDims <- ncol(points)
    if (nDims != 2 && nDims != 3)
        report(OL$Error, "Points must be two or three dimensional")
    
    result <- .Call("cp_transform_R", .fixTypes(controlPointImage), .fixTypes(target), points, as.logical(nearest), PACKAGE="RNiftyReg")
    
    newPoints <- sapply(seq_len(nrow(points)), function(i) {
        if (length(result[[i]]) == nDims)
            return (result[[i]])
        else
        {
            data <- as.data.frame(matrix(result[[i]], ncol=2*nDims, byrow=TRUE))
            if (nDims == 2)
            {
                colnames(data) <- c("sx", "sy", "tx", "ty")
                fit <- lm(cbind(tx,ty) ~ bs(sx) * bs(sy), data=data)
                return (drop(predict(fit, data.frame(sx=points[i,1],sy=points[i,2]))))
            }
            else
            {
                colnames(data) <- c("sx", "sy", "sz", "tx", "ty", "tz")
                fit <- lm(cbind(tx,ty,tz) ~ bs(sx) * bs(sy) * bs(sz), data=data)
                return (drop(predict(fit, data.frame(sx=points[i,1],sy=points[i,2],sz=points[i,3]))))
            }
        }
    })
    
    dimnames(newPoints) <- NULL
    newPoints <- drop(t(newPoints))
    
    return (newPoints)
}
