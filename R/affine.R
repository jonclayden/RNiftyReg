readAffine <- function (fileName, type = NULL)
{
    if (!is.null(type))
        type <- match.arg(tolower(type), c("niftyreg","fsl"))
    
    lines <- readLines(fileName)
    typeLine <- (lines %~% "\\# affineType\\: \\w+")
    if (!any(typeLine) && is.null(type))
        report(OL$Error, "Type is not stored in the file - it must be specified")
    else if (is.null(type))
        type <- match.arg(tolower(sub("\\# affineType\\: (\\w+)", "\\1", lines[typeLine][1], perl=TRUE)), c("niftyreg","fsl"))
    
    connection <- textConnection(lines[!typeLine])
    affine <- as.matrix(read.table(connection))
    close(connection)
    
    if (!isTRUE(all.equal(dim(affine), c(4,4))))
        report(OL$Error, "The specified file does not contain a 4x4 affine matrix")
    
    attr(affine, "affineType") <- type
    return (affine)
}

writeAffine <- function (affine, fileName)
{
    if (!is.matrix(affine) || !isTRUE(all.equal(dim(affine), c(4,4))))
        report(OL$Error, "Specified affine matrix is not valid")
    
    lines <- apply(format(affine,scientific=FALSE),1,implode,sep="  ")
    lines <- c(paste("# affineType:",attr(affine,"affineType"),sep=" "), lines)
    writeLines(lines, fileName)
}

convertAffine <- function (affine, source, target, newType = c("niftyreg","fsl"), currentType = NULL)
{
    if (!is.matrix(affine) || !isTRUE(all.equal(dim(affine), c(4,4))))
        report(OL$Error, "Specified affine matrix is not valid")
    if (!is.nifti(source) || !is.nifti(target))
        report(OL$Error, "Source and target images must be specified as \"nifti\" objects")
    
    newType <- match.arg(newType)
    
    if (is.null(currentType))
    {
        currentType <- attr(affine, "affineType")
        if (is.null(currentType))
            report(OL$Error, "The current affine type was not specified and is not stored with the matrix")
    }
    else
        currentType <- match.arg(currentType, c("niftyreg","fsl"))
    
    if (newType == currentType)
        return (affine)
    else
    {
        sourceXform <- xformToAffine(source, useQuaternionFirst=FALSE)
        targetXform <- xformToAffine(target, useQuaternionFirst=FALSE)
        sourceScaling <- diag(c(sqrt(colSums(sourceXform[1:3,1:3]^2)), 1))
        targetScaling <- diag(c(sqrt(colSums(targetXform[1:3,1:3]^2)), 1))
        
        if (newType == "fsl")
            newAffine <- targetScaling %*% solve(targetXform) %*% solve(affine) %*% sourceXform %*% solve(sourceScaling)
        else
            newAffine <- sourceXform %*% solve(sourceScaling) %*% solve(affine) %*% targetScaling %*% solve(targetXform)
        
        attr(newAffine, "affineType") <- newType
        return (newAffine)
    }
}

decomposeAffine <- function (affine, source = NULL, target = NULL, type = NULL)
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
    
    # Full matrix is rotationX %*% rotationY %*% rotationZ %*% skew %*% scale
    submatrix <- affine[1:3,1:3]
    sm <- list(x=submatrix[,1], y=submatrix[,2], z=submatrix[,3])
    xLength <- sqrt(sum(sm$x^2))
    yLength <- sqrt((sm$y %*% sm$y) - (sm$x %*% sm$y)^2 / xLength^2)
    xyProj <- (sm$x %*% sm$y) / (xLength * yLength)
    xNorm <- sm$x / xLength
    yNorm <- (sm$y / yLength) - (xyProj * xNorm)
    zLength <- sqrt((sm$z %*% sm$z) - (xNorm %*% sm$z)^2 - (yNorm %*% sm$z)^2)
    xzProj <- (xNorm %*% sm$z) / zLength
    yzProj <- (yNorm %*% sm$z) / zLength
    
    scales <- c(xLength, yLength, zLength)
    scaleMatrix <- diag(scales)
    skews <- c(xyProj, xzProj, yzProj)
    skewMatrix <- diag(3)
    skewMatrix[c(4,7,8)] <- skews
    translation <- affine[1:3,4]
    
    rotationMatrix <- submatrix %*% solve(scaleMatrix) %*% solve(skewMatrix)
    pitchAngle <- asin(-rotationMatrix[1,3])
    if (cos(pitchAngle) < 1e-4)
    {
        # Degenerate case (Gimbal lock) - fix yaw angle at zero
        rollAngle <- atan2(-rotationMatrix[3,2], rotationMatrix[2,2])
        yawAngle <- 0
    }
    else
    {
        rollAngle <- atan2(rotationMatrix[2,3], rotationMatrix[3,3])
        yawAngle <- atan2(rotationMatrix[1,2], rotationMatrix[1,1])
    }
    angles <- c(rollAngle, pitchAngle, yawAngle)
    
    names(translation) <- letters[24:26]
    names(scales) <- letters[24:26]
    names(skews) <- c("xy", "xz", "yz")
    names(angles) <- c("roll", "pitch", "yaw")
    
    return (list(scaleMatrix=scaleMatrix, skewMatrix=skewMatrix, rotationMatrix=rotationMatrix, translation=translation, scales=scales, skews=skews, angles=angles))
}
