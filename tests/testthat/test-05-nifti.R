context("Reading and writing NIfTI files")

test_that("NIfTI files can be read and written", {
    imagePath <- system.file("extdata", "epi_t2.nii.gz", package="RNiftyReg")
    tempPath <- tempfile()
    
    expect_that(dim(readNifti(imagePath,internal=FALSE)), equals(c(96L,96L,60L)))
    expect_that(dim(readNifti(imagePath,internal=TRUE)), equals(c(96L,96L,60L)))
    
    image <- readNifti(imagePath)
    expect_that(pixunits(image), equals(c("mm","s")))
    writeNifti(image, tempPath)
    expect_that(pixdim(readNifti(tempPath)), equals(c(2.5,2.5,2.5)))
    unlink(tempPath)
})
