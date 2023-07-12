options(RNiftyReg.threads=2L)

t2 <- readNifti(system.file("extdata","epi_t2.nii.gz",package="RNiftyReg"))
t1 <- readNifti(system.file("extdata","flash_t1.nii.gz",package="RNiftyReg"))
mni <- readNifti(system.file("extdata","mni_brain.nii.gz",package="RNiftyReg"))

t2_to_t1 <- readAffine(system.file("extdata","affine.txt",package="RNiftyReg"), t2, t1)
t1_to_mni <- readNifti(system.file("extdata","control.nii.gz",package="RNiftyReg"), t1, mni)

deformation <- deformationField(t2_to_t1, jacobian=TRUE)
expect_equal(round(worldToVoxel(as.array(deformation)[34,49,64,1,], t2)), c(40,40,20))
expect_equal(as.array(jacobian(deformation))[34,49,64], prod(diag(t2_to_t1)), tolerance=0.05)

expect_equal(applyTransform(t2_to_t1,c(40,40,20),nearest=TRUE), c(34,49,64))
expect_equal(class(applyTransform(t2_to_t1,t2,internal=TRUE))[1], "internalImage")

rdsFile <- tempfile(fileext="rds")
saveTransform(t2_to_t1, rdsFile)
reloadedTransform <- loadTransform(rdsFile)
expect_equal(applyTransform(reloadedTransform,c(40,40,20),nearest=TRUE), c(34,49,64))
expect_equivalent(applyTransform(t2_to_t1,t2), applyTransform(reloadedTransform,t2))

if (tolower(Sys.info()[["sysname"]]) != "sunos") {
    point <- applyTransform(t2_to_t1, c(40,40,20), nearest=FALSE)
    expect_equal(applyTransform(t1_to_mni,point,nearest=TRUE), c(33,49,24))
    expect_equal(round(applyTransform(t1_to_mni,point,nearest=FALSE)), c(33,49,24))
    
    saveTransform(t1_to_mni, rdsFile)
    reloadedTransform <- loadTransform(rdsFile)
    expect_equal(applyTransform(reloadedTransform,point,nearest=TRUE), c(33,49,24))
    
    # Different z-value due to double-rounding
    expect_equal(applyTransform(t1_to_mni,t1,interpolation=0)[33,49,25], t1[34,49,64])
    
    # Extract affine embedded in extensions
    expect_inherits(asAffine(t1_to_mni), "affine")
    
    t2_to_mni <- composeTransforms(t2_to_t1, t1_to_mni)
    expect_equal(applyTransform(t2_to_mni,c(40,40,20),nearest=TRUE), c(33,49,24))
    
    t2_to_t1_half <- halfTransform(t2_to_t1)
    expect_equivalent(composeTransforms(t2_to_t1_half,t2_to_t1_half), t2_to_t1)
    
    # Use an identity transform to ensure that the target is right
    t1_to_mni_half <- halfTransform(t1_to_mni)
    mniIdentity <- buildAffine(source=mni)
    t1_to_mni_reconstructed <- composeTransforms(t1_to_mni_half, t1_to_mni_half, mniIdentity)
    expect_equal(applyTransform(t1_to_mni_reconstructed,point,nearest=TRUE), c(33,49,24))
}
