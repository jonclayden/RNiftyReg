context("Affine matrix operations")

test_that("Affine operations work", {
    source <- readNifti(system.file("extdata","epi_t2.nii.gz",package="RNiftyReg"), internal=FALSE)
    target <- readNifti(system.file("extdata","flash_t1.nii.gz",package="RNiftyReg"), internal=FALSE)
    
    affine <- readAffine(system.file("extdata","affine.txt",package="RNiftyReg"), source, target)
    fslAffine <- readAffine("flirt.mat", source, target, "fsl")
    
    # FSL and NiftyReg transforms are fairly similar for the same source and target images
    expect_that(fslAffine, equals(affine,tolerance=0.1,check.attributes=FALSE))
    expect_that(isAffine(affine), equals(TRUE))
    expect_that(print(affine), prints_text("origin",fixed=TRUE))
    expect_that(invertAffine(invertAffine(affine)), is_equivalent_to(affine))
    expect_that(buildAffine(decomposeAffine(affine),source=source,target=target), is_equivalent_to(affine))
    
    expect_that(buildAffine(angles=c(0,0,pi/4),source=source)[,4], equals(c(0,0,0,1)))
    expect_that(round(buildAffine(angles=c(0,0,pi/4),source=source,anchor="centre")[,4]), equals(c(18,5,0,1)))
    
    origin <- function(image) worldToVoxel(c(0,0,0), image)
    xfm <- buildAffine(scales=c(2,2,2), source=source)
    expect_that(origin(attr(xfm,"target")), equals((origin(attr(xfm,"source"))-1)*2+1))
})
