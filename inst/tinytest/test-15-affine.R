source <- readNifti(system.file("extdata","epi_t2.nii.gz",package="RNiftyReg"), internal=FALSE)
target <- readNifti(system.file("extdata","flash_t1.nii.gz",package="RNiftyReg"), internal=FALSE)

affine <- readAffine(system.file("extdata","affine.txt",package="RNiftyReg"), source, target)
fslAffine <- readAffine("flirt.mat", source, target, "fsl")

# FSL and NiftyReg transforms are fairly similar for the same source and target images
expect_equivalent(fslAffine, affine, tolerance=0.1)
expect_equal(isAffine(affine), TRUE)
expect_stdout(print(affine), "origin")
expect_equivalent(invertAffine(invertAffine(affine)), affine)
expect_equivalent(buildAffine(decomposeAffine(affine),source=source,target=target), affine)

expect_equal(buildAffine(angles=c(0,0,pi/4),source=source)[,4], c(0,0,0,1))
expect_equal(round(buildAffine(angles=c(0,0,pi/4),source=source,anchor="centre")[,4]), c(18,5,0,1))

xfm <- buildAffine(scales=c(2,2,2), source=source)
expect_equal(RNifti::origin(attr(xfm,"target")), (RNifti::origin(attr(xfm,"source"))-1)*2+1)

affineFile <- tempfile(fileext="txt")
writeAffine(affine, affineFile)
reloadedAffine <- readAffine(affineFile, source, target)
expect_inherits(reloadedAffine, "affine")
expect_equivalent(affine, reloadedAffine)
