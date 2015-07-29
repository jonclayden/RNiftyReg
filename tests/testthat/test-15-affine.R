context("Affine matrix operations")

test_that("Affine operations work", {
    source <- readNifti(system.file("extdata","epi_t2.nii.gz",package="RNiftyReg"), internal=FALSE)
    target <- readNifti(system.file("extdata","flash_t1.nii.gz",package="RNiftyReg"), internal=FALSE)
    
    affine <- readAffine(system.file("extdata","affine.txt",package="RNiftyReg"), source, target)
    fslAffine <- readAffine("flirt.mat", source, target, "fsl")
    
    # FSL and NiftyReg transforms are fairly similar for the same source and target images
    expect_that(fslAffine, equals(affine,tolerance=0.1,check.attributes=FALSE))
    expect_that(isAffine(affine), equals(TRUE))
    expect_that(invertAffine(invertAffine(affine)), is_equivalent_to(affine))
    expect_that(buildAffine(decomposeAffine(affine),source=source,target=target), is_equivalent_to(affine))
})
