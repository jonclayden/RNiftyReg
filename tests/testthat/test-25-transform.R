context("Applying transformations")

test_that("Existing transformations can be applied and combined", {
    t2 <- readNifti(system.file("extdata","epi_t2.nii.gz",package="RNiftyReg"))
    t1 <- readNifti(system.file("extdata","flash_t1.nii.gz",package="RNiftyReg"))
    mni <- readNifti(system.file("extdata","mni_brain.nii.gz",package="RNiftyReg"))
    
    t2_to_t1 <- readAffine(system.file("extdata","affine.txt",package="RNiftyReg"), t2, t1)
    t1_to_mni <- readNifti(system.file("extdata","control.nii.gz",package="RNiftyReg"), t1, mni)
    
    deformation <- deformationField(t2_to_t1, jacobian=TRUE)
    expect_that(round(worldToVoxel(as.array(deformation)[34,49,64,1,], t2)), equals(c(40,40,20)))
    expect_that(as.array(jacobian(deformation))[34,49,64], equals(prod(diag(t2_to_t1)),tolerance=0.05))
    
    expect_that(applyTransform(t2_to_t1,c(40,40,20),nearest=TRUE), equals(c(34,49,64)))
    expect_that(class(applyTransform(t2_to_t1,t2,internal=TRUE))[1], equals("internalImage"))
    
    skip_on_os("solaris")
    
    point <- applyTransform(t2_to_t1, c(40,40,20), nearest=FALSE)
    expect_that(applyTransform(t1_to_mni,point,nearest=TRUE), equals(c(33,49,24)))
    expect_that(round(applyTransform(t1_to_mni,point,nearest=FALSE)), equals(c(33,49,24)))
    
    # Different z-value due to double-rounding
    expect_that(applyTransform(t1_to_mni,t1,interpolation=0)[33,49,25], equals(t1[34,49,64]))
    
    t2_to_mni <- composeTransforms(t2_to_t1, t1_to_mni)
    expect_that(applyTransform(t2_to_mni,c(40,40,20),nearest=TRUE), equals(c(33,49,24)))
    
    t2_to_t1_half <- halfTransform(t2_to_t1)
    expect_that(composeTransforms(t2_to_t1_half,t2_to_t1_half), is_equivalent_to(t2_to_t1))
    
    t1_to_mni_half <- halfTransform(t1_to_mni)
    t1_to_mni_reconstructed <- composeTransforms(t1_to_mni_half, t1_to_mni_half)
    expect_that(applyTransform(t1_to_mni_reconstructed,point,nearest=TRUE), equals(c(33,49,24)))
})
