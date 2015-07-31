context("Applying transformations")

test_that("Existing transformations can be applied and combined", {
    t2 <- readNifti(system.file("extdata","epi_t2.nii.gz",package="RNiftyReg"))
    t1 <- readNifti(system.file("extdata","flash_t1.nii.gz",package="RNiftyReg"))
    mni <- readNifti(system.file("extdata","mni_brain.nii.gz",package="RNiftyReg"))
    
    t2_to_t1 <- readAffine(system.file("extdata","affine.txt",package="RNiftyReg"), t2, t1)
    t1_to_mni <- readNifti(system.file("extdata","control.nii.gz",package="RNiftyReg"), t1, mni)
    
    expect_that(applyTransform(t2_to_t1,c(40,40,20),nearest=TRUE), equals(c(34,49,64)))
    expect_that(applyTransform(t1_to_mni,c(34,49,64),nearest=TRUE), equals(c(33,49,25)))
    
    # Result here is slightly different, due to rounding error in the intermediate location
    t2_to_mni <- composeTransforms(t2_to_t1, t1_to_mni)
    expect_that(applyTransform(t2_to_mni,c(40,40,20),nearest=TRUE), equals(c(33,49,24)))
    
    t2_to_t1_half <- halfTransform(t2_to_t1)
    expect_that(composeTransforms(t2_to_t1_half,t2_to_t1_half), is_equivalent_to(t2_to_t1))
    
    t1_to_mni_half <- halfTransform(t1_to_mni)
    t1_to_mni_reconstructed <- composeTransforms(t1_to_mni_half, t1_to_mni_half)
    expect_that(applyTransform(t1_to_mni_reconstructed,c(34,49,64),nearest=TRUE), equals(c(33,49,25)))
})
