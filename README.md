# RNiftyReg: Nifty Registration in R

The `RNiftyReg` package is an R-native interface to the [NiftyReg image registration library](http://sourceforge.net/projects/niftyreg/) developed within the Translational Imaging Group at University College London. The package incorporates the library, so it does not need to be installed separately, and it replaces the NiftyReg command-line front-end with a direct, in-memory bridge to R, based on [Rcpp](http://www.rcpp.org).

This `README` file primarily covers version 2.0.0 of the package and later. The interface was substantially reworked in that version to make it more natural and less verbose, and earlier versions are incompatible. Information on moving from prior versions of `RNiftyReg` to 2.x is included at the end of this file.

## Contents

- [Reading and writing images](#reading-and-writing-images)
- [Image registration](#image-registration)
- [Applying transformations](#applying-transformations)
- [RNiftyReg internals](#rniftyreg-internals)
- [Upgrading to RNiftyReg 2.x](#upgrading-to-rniftyreg-2x)

## Reading and writing images

`RNiftyReg` may be used to register and manipulate two and three dimensional images of any sort, although its origins are in medical imaging. Medical images in the standard [NIfTI-1 format](http://nifti.nimh.nih.gov/nifti-1) may be read into R using the `readNifti` function.

```r
image <- readNifti(system.file("extdata", "epi_t2.nii.gz", package="RNiftyReg"))
```

This image is an R array with some additional attributes containing information such as its dimensions and the size of its pixels (or voxels, in this case, since it is a 3D image). There are auxiliary functions for extracting this information: the standard `dim()`, plus `pixdim()` and `pixunits()`.

```r
dim(image)
# [1] 96 96 60

pixdim(image)
# [1] 2.5 2.5 2.5

pixunits(image)
# [1] "mm" "s"
```

So this image is of size 96 x 96 x 60 voxels, with each voxel representing 2.5 x 2.5 x 2.5 mm in real space. (The temporal unit, seconds here, only applies to the fourth dimension, if it is present.) An image can be written back to NIfTI-1 format using the complementary `writeNifti` function.

```r
writeNifti(image, "file.nii.gz")
```

As mentioned above, however, images do not have to be in NIfTI-1 format. Any numeric matrix or array can be used, and standard image formats such as JPEG and PNG can be read in using additional packages. For example,

```r
library(png)
image <- readPNG(system.file("extdata", "house.png", package="RNiftyReg"))
```

The utility functions mentioned above can still be applied, but defaults are returned where necessary.

```r
dim(image)
# [1] 182 261

pixdim(image)
# [1] 1 1

pixunits(image)
# [1] "Unknown"
```

## Image registration

## Applying transformations

## RNiftyReg internals

## Upgrading to RNiftyReg 2.x
