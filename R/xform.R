xformToAffine <- function (image, useQuaternionFirst = TRUE)
{
    if (!require("oro.nifti"))
        report(OL$Error, "The \"oro.nifti\" package is required")
    if (!is.nifti(image))
        report(OL$Error, "The specified image is not a \"nifti\" object")
    
    # With no information, assume Analyze orientation and zero origin
    if (image@qform_code <= 0 && image@sform_code <= 0)
        return (diag(c(-1, 1, 1, 1)))
    else if ((useQuaternionFirst && image@qform_code > 0) || image@sform_code <= 0)
    {
        matrix <- diag(4)
        matrix[1:3,4] <- c(image@qoffset_x, image@qoffset_y, image@qoffset_z)
        q <- c(image@quatern_b, image@quatern_c, image@quatern_d)
        
        if (sum(q^2) == 1)
            q <- c(0, q)
        else
            q <- c(sqrt(1 - sum(q^2)), q)
        
        matrix[1:3,1:3] <- c(q[1]*q[1] + q[2]*q[2] - q[3]*q[3] - q[4]*q[4],
                             2*q[2]*q[3] + 2*q[1]*q[4],
                             2*q[2]*q[4] - 2*q[1]*q[3],
                             2*q[2]*q[3] - 2*q[1]*q[4],
                             q[1]*q[1] + q[3]*q[3] - q[2]*q[2] - q[4]*q[4],
                             2*q[3]*q[4] + 2*q[1]*q[2],
                             2*q[2]*q[4] + 2*q[1]*q[3],
                             2*q[3]*q[4] - 2*q[1]*q[2],
                             q[1]*q[1] + q[4]*q[4] - q[3]*q[3] - q[2]*q[2])
        
        # The qfactor should be stored as 1 or -1, but the NIfTI standard says
        # 0 should be treated as 1
        # This formulation does that (the 0.1 is arbitrary)
        qfactor <- sign(image@pixdim[1] + 0.1)
        matrix[1:3,1:3] <- matrix[1:3,1:3] * rep(c(abs(image@pixdim[2:3]), qfactor*abs(image@pixdim[4])), each=3)
        
        return (matrix)
    }
    else
        return (rbind(image@srow_x, image@srow_y, image@srow_z, c(0,0,0,1)))
}
