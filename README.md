# RNiftyReg: Nifty Registration in R

The `RNiftyReg` package is an R-native interface to the [NiftyReg image registration library](http://sourceforge.net/projects/niftyreg/) developed within the Translational Imaging Group at University College London. The package incorporates the library, so it does not need to be installed separately, and it replaces the NiftyReg command-line front-end with a direct, in-memory bridge to R, based on [Rcpp](http://www.rcpp.org).

This `README` file primarily covers version 2.0.0 of the package and later. The interface was substantially reworked in that version to make it more natural and less verbose, and earlier versions are incompatible. Information on moving from prior versions of `RNiftyReg` to 2.x is included at the end of this file.

The package can be installed from CRAN, or the latest development version obtained directly from GitHub:

```r
## install.packages("devtools")
devtools::install_github("jonclayden/RNiftyReg")
```

The [`mmand` package](https://github.com/jonclayden/mmand) for image processing may also be useful, and is used in some of the examples below.

## Contents

- [Reading and writing images](#reading-and-writing-images)
- [Image registration](#image-registration)
- [Applying transformations](#applying-transformations)
- [RNiftyReg internals](#rniftyreg-internals)
- [Upgrading to RNiftyReg 2.x](#upgrading-to-rniftyreg-2x)

## Reading and writing images

`RNiftyReg` may be used to register and manipulate two and three dimensional images of any sort, although its origins are in medical imaging. Medical images in the standard [NIfTI-1 format](http://nifti.nimh.nih.gov/nifti-1) may be read into R using the `readNifti` function.

```r
library(RNiftyReg)
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
## install.packages("jpeg")
library(jpeg)
image <- readJPEG(system.file("extdata", "house_colour_large.jpg", package="RNiftyReg"))
```

The utility functions mentioned above can still be applied, but defaults are returned where necessary.

```r
dim(image)
# [1] 363 523   3

pixdim(image)
# [1] 1 1 1

pixunits(image)
# [1] "Unknown"
```

The pixel dimensions can be set manually if required.

```r
pixdim(image) <- c(2, 2)
```

## Image registration

Once two or more images have been read into R, they can be registered. Registration is the dual operation of searching a space of transformations for the best way to align two images, and then resampling one image onto the grid of the other. The images need not be the same size or in the same orientation.

There are two main classes of transformation available: linear and nonlinear. Linear, and specifically *affine*, transforms can represent translation, scaling and rotation in 2D or 3D space. They have up to 12 degrees of freedom and are appropriate to capture global shifts between images. Nonlinear transformations have many more degrees of freedom, and can capture localised distortions between images, but they are more time-consuming to estimate and more complex to work with.

Some sample 3D medical images are included with the package. We begin by registering two brain scan images, of the same person, with different contrasts. First we read them in, and then we pass them to the package's core registration function, `niftyreg`.

```r
source <- readNifti(system.file("extdata", "epi_t2.nii.gz", package="RNiftyReg"))
target <- readNifti(system.file("extdata", "flash_t1.nii.gz", package="RNiftyReg"))

result <- niftyreg(source, target)
```

The last command will take a few seconds to complete. The `result` is a list with a number of components, the most important of which is the resampled source image in target space, `result$image`.

By default the transformation will be an affine matrix. If we want to allow for a nonlinear transformation, with scope for local deformations between one space and the other, we can perform an additional registration as follows.

```r
result <- niftyreg(source, target, scope="nonlinear", init=forward(result))
```

Notice the `scope` argument, and also the fact that we use the result of the previous linear registration to initialise the nonlinear one. (The `forward` function extracts the forward transformation from the previous registration.) This should result in reduced convergence time, since the affine transformation provides a first approximation for the nonlinear registration algorithm. Nevertheless, this algorithm will generally take longer to complete.

## Applying transformations

Once a transformation between two images has been established through registration, it can be extracted, applied to another image or pixel coordinates, or manipulated. Transformations can also be read from or written to file, or created from scratch. Registration is by default symmetric, meaning that forward and reverse transformations are calculated simultaneously. These can be extracted using the `forward` and `reverse` functions.

Let's use a simple image by way of example. We will need the `jpeg` package to read it in, and the `mmand` package (version 1.2.0 or later) to visualise it.

```r
## install.packages(c("jpeg","mmand"))
library(jpeg); library(mmand)

house <- readJPEG(system.file("extdata", "house_colour_large.jpg", package="RNiftyReg"))
display(house)
```

![House photo](README-house.jpg)

Clearly this is a colour image, with red, green and blue channels, but NiftyReg cannot work directly with it in this format, so we convert it to greyscale by simply averaging across channels:

```r
house_bw <- apply(house, 1:2, mean)
display(house_bw)
```

![Greyscale house photo](README-house-bw.jpg)

Now, instead of registering the image to another image, let's create a simple affine transformation that applies a skew to the image.

```r
affine <- buildAffine(skews=0.1, source=house_bw, target=house_bw)
print(affine)
# NiftyReg affine matrix:
#  1.0  -0.1   0.0   0.0
#  0.0   1.0   0.0   0.0
#  0.0   0.0   1.0   0.0
#  0.0   0.0   0.0   1.0
# Source origin: (1, 1, 1)
# Target origin: (1, 1, 1)
```

So, this is a diagonal matrix with just a single off-diagonal element, which produces the skew effect. Let's apply it to the image using the important `applyTransform` function, and see the effect.

```r
house_skewed_bw <- applyTransform(affine, house_bw)
display(house_skewed_bw)
```

![Skewed greyscale house photo](README-house-skewed-bw.jpg)

Notice that we can also transform the colour image; in this case, `applyTransform` applies the transform to each channel individually to produce the final result:

```r
house_skewed <- applyTransform(affine, house)
display(house_skewed)
```

![Skewed colour house photo](README-house-skewed.jpg)

Moreover, we can transform a pixel coordinate into the space of the skewed image:

```r
applyTransform(affine, c(182,262))
# [1] 208.1 262.0
```

Notice that the skew changes the first coordinate (in the up-down direction), but not the second (in the left-right direction).

Finally, we can register the original image to the skewed one, to recover the transformation:

```r
result <- niftyreg(house_bw, house_skewed_bw, scope="affine")
print(forward(result))
# NiftyReg affine matrix:
#  1.000845194  -0.099988878   0.000000000  -0.279240519
# -0.003345727   1.000585556   0.000000000   0.543151319
#  0.000000000   0.000000000   1.000000000   0.000000000
#  0.000000000   0.000000000   0.000000000   1.000000000
# Source origin: (1, 1, 1)
# Target origin: (1, 1, 1)
```

Notice that the estimated transformation closely approximates the generative one, with the element in row 1, column 2 being very close to -0.1. We can decompose this estimated transformation and recover the skew component:

```r
decomposeAffine(forward(result))$skews
#        xy        xz        yz 
# 0.1032827 0.0000000 0.0000000 
```

Two other manipulations can be helpful to know about. The first is calculating a half-transform, which can be used to transform the image into a space halfway to the target. For example, using our registration result from above,

```r
half_xfm <- halfTransform(forward(result))
display(applyTransform(half_xfm, house))
```

![Half-skewed colour house photo](README-house-halfskewed.jpg)

This results in half of the skew effect being applied. Finally, the `composeTransforms` function allows the effects of two transforms to be combined together. Combining a half-transform with itself will result in the original full transform.

```r
all.equal(forward(result), composeTransforms(half_xfm,half_xfm), check.attributes=FALSE)
# TRUE
```

## RNiftyReg internals

**NB**: This section is technical and will not be of interest to a lot of users. Please feel free to skip it!

NiftyReg has its origins in medical imaging, and as such it natively uses a data structure specific to the [NIfTI-1 file format](http://nifti.nimh.nih.gov/nifti-1). `RNiftyReg` uses the standard NIfTI-1 library to read, write and manipulate these data structures, which is fast and well-tested.

Importantly, the package does not fully duplicate the NIfTI-1 structure's contents in R-visible objects. Instead, it passes key metadata back to R, such as the image dimensions and pixel dimensions, and it also passes back the pixel values where they are needed. It also passes back an [external pointer](http://r-manuals.flakery.org/R-exts.html#External-pointers-and-weak-references) to the native data structure, which is stored in an attribute. This pointer is dereferenced whenever the object is passed back to the C++ code, thereby avoiding unnecessary duplication and ensuring that all metadata remains intact. The full NIfTI-1 header can be obtained using the `dumpNifti` R function, if it is needed.

This arrangement is efficient and generally works well, but there is the possibility that the image data could get out of sync between the C++ data structure and the R array. Many R operations strip attributes—in which case the external pointer will be removed and the internal structure will be built again when necessary—but some do not. In those cases, the `updateNifti` function should be called on the object to synchronise the internal version of the data with R-visible array. This happens automatically in the registration functions, and when writing a NIfTI-1 image out with `writeNifti`, but in other places it may be necessary to call the function explicitly. If the object has no external pointer attribute, or if the metadata needs to be overwritten, a template image can be supplied.

## Upgrading to RNiftyReg 2.x

`RNiftyReg` 2.0.0 is a more-or-less complete rewrite of the package, with the goals of simplifying both the package's dependencies and its usage. The upstream NiftyReg code has also been updated. The core changes are

- The `oro.nifti` package is no longer needed, nor used for reading and writing NIfTI files (`RNiftyReg` now offers `readNifti` and `writeNifti`, which are much faster). However, objects of S4 class `nifti` can still be used with the package if desired. Functions return either plain R arrays with attributes or bare-bones `internalImage` objects, which contain only some basic metadata and a pointer to a C-level data structure.
- There are new functions for halving a transform (`halfTransform`), composing two transforms (`composeTransforms`), and building an affine transform from scratch (`buildAffine`).
- Registration is now symmetric by default (for both linear and nonlinear), a newer symmetric nonlinear approach is now used, and default cost function weights have been tweaked. Therefore, the arguments to the core `niftyreg` function, and its linear and nonlinear special cases, have changed in both name and defaults. See `?niftyreg` and related help pages for details.
- It is no longer necessary to use functions specific to transform type to perform many operations. For example, the work of the old `applyAffine`, `applyControlPoints`, `transformWithAffine` and `transformWithControlPoints` functions is done by the flexible new `applyTransform` function. The forward and reverse transforms can always be obtained from a registration using the new `forward` and `reverse` functions, no matter what their type is. However, some affine-only functions, such as `decomposeAffine`, retain their names.
- The `affineType` attribute has gone, and `convertAffine` is no longer a user-visible function. All affine matrices are stored using the NiftyReg convention. FSL-FLIRT affines can still be read in, but they are converted to NiftyReg convention immediately. In addition, source and target image information is attached to the transforms in attributes, and so does not need to be specified in most function calls.

However, it should still be possible to read and use transformations created using `RNiftyReg` 1.x.
