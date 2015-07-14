applyTransform <- function (transform, x, interpolation = 3L)
{
    source <- attr(transform, "source")
    target <- attr(transform, "target")
    
    if (isAffine(transform, strict=TRUE))
    {
        # The argument looks like a suitable image
        if (isImage(x,TRUE) && isTRUE(all.equal(dim(x),dim(source))))
            return (niftyreg.linear(x, target, "affine", init=transform, nLevels=0L, interpolation=interpolation, verbose=FALSE, estimateOnly=FALSE))
        else if ((is.matrix(x) && ncol(x) == length(dim(source))) || length(x) == length(dim(source)))
        {
            fslAffine <- convertAffine(transform, source, target, "fsl")
            points <- voxelToWorld(x, source, simple=TRUE)
            newPoints <- applyAffine(points, affine)
            newPoints <- worldToVoxel(newPoints, target, simple=TRUE)
            return (newPoints)
        }
        else
            report(OL$Error, "Object to transform should be a suitable image or matrix of points")
    }
    else if (isImage(transform, FALSE))
    {
        if (isImage(x,TRUE) && isTRUE(all.equal(dim(x),dim(source))))
            return (niftyreg.nonlinear(x, target, init=transform, nLevels=0L, interpolation=interpolation, verbose=FALSE, estimateOnly=FALSE))
        else if ((is.matrix(x) && ncol(x) == length(dim(source))) || length(x) == length(dim(source)))
        {
            points <- voxelToWorld(points, source)
            
            if (!is.matrix(points))
                points <- matrix(points, nrow=1)
            
            nDims <- ncol(points)
            if (nDims != 2 && nDims != 3)
                report(OL$Error, "Points must be two or three dimensional")
            
            result <- .Call("cp_transform_R", .fixTypes(transform), .fixTypes(target), points, as.logical(nearest), PACKAGE="RNiftyReg")
            
            newPoints <- sapply(seq_len(nrow(points)), function(i) {
                if (length(result[[i]]) == nDims)
                    return (result[[i]])
                else
                {
                    data <- as.data.frame(matrix(result[[i]], ncol=2*nDims, byrow=TRUE))
                    if (nDims == 2)
                    {
                        colnames(data) <- c("sx", "sy", "tx", "ty")
                        fit <- lm(cbind(tx,ty) ~ splines::bs(sx) * splines::bs(sy), data=data)
                        return (drop(predict(fit, data.frame(sx=points[i,1],sy=points[i,2]))))
                    }
                    else
                    {
                        colnames(data) <- c("sx", "sy", "sz", "tx", "ty", "tz")
                        fit <- lm(cbind(tx,ty,tz) ~ splines::bs(sx) * splines::bs(sy) * splines::bs(sz), data=data)
                        return (drop(predict(fit, data.frame(sx=points[i,1],sy=points[i,2],sz=points[i,3]))))
                    }
                }
            })
            
            dimnames(newPoints) <- NULL
            newPoints <- drop(t(newPoints))
            return (newPoints)
        }
        else
            report(OL$Error, "Object to transform should be a suitable image or matrix of points")
    }
    else
        report(OL$Error, "Specified transform is not valid")
}
