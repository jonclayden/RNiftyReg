context("Reading and writing NIfTI files")

test_that("NIfTI files can be read and written", {
    sourcePath <- system.file("extdata", "source.nii.gz", package="RNiftyReg")
    tempPath <- tempfile()
    
    expect_that(dim(readNifti(sourcePath,internal=FALSE)), equals(c(96L,96L,25L)))
    expect_that(dim(readNifti(sourcePath,internal=TRUE)), equals(c(96L,96L,25L)))
    
    image <- readNifti(sourcePath)
    expect_that(pixunits(image), equals(c("mm","s")))
    writeNifti(image, tempPath)
    expect_that(pixdim(readNifti(tempPath)), equals(c(2.5,2.5,5)))
    unlink(tempPath)
})
