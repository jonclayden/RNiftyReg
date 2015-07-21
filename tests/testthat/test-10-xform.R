context("NIfTI sform/qform operations")

test_that("NIfTI sform/qform operations work", {
    image <- readNifti(system.file("extdata", "source.nii.gz", package="RNiftyReg"))
    
    expect_that(diag(xform(image)), equals(c(-2.5,2.5,5,1)))
    expect_that(diag(xform(image,useQuaternionFirst=FALSE)), equals(c(-2.5,2.5,5,1)))
    
    point <- c(40, 40, 10)
    expect_that(voxelToWorld(point,image), equals(c(-97.5,97.5,45.0)))
    expect_that(voxelToWorld(point,image,simple=TRUE), equals(c(97.5,97.5,45.0)))
    expect_that(worldToVoxel(voxelToWorld(point,image),image), equals(point))
    expect_that(worldToVoxel(voxelToWorld(point,image,simple=TRUE),image,simple=TRUE), equals(point))
})
