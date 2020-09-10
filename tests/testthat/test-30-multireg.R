context("Multiple registration")

test_that("Multiple registration works", {
    options(RNiftyReg.threads=2L)

    if (system.file(package="loder") == "")
        skip("The \"loder\" package is not available")
    else
    {
        house <- loder::readPng(system.file("extdata","house.png",package="RNiftyReg"))
        affine <- buildAffine(skews=0.1, source=house, target=house)
        skewedHouse <- applyTransform(affine, house)
        
        colourHouse <- loder::readPng(system.file("extdata","house_colour.png",package="RNiftyReg"))
        skewedColourHouse <- applyTransform(affine, colourHouse)
        
        # Rec. 709 luma RGB-to-greyscale coefficients (as used by ImageMagick)
        expect_equal(skewedHouse[66,76], weighted.mean(skewedColourHouse[66,76,],c(0.2126,0.7152,0.0722)), tolerance=0.05)
        
        skip_on_os("solaris")
        
        reg <- niftyreg(skewedColourHouse, house, symmetric=FALSE)
        expect_equal(forward(reg)[1,2], 0.1, tolerance=0.05)
        
        skip_on_cran()
        skip_on_travis()
        
        reg <- niftyreg(skewedColourHouse, house, scope="nonlinear", init=forward(reg), symmetric=FALSE, maxIterations=300)
        expect_equal(dim(forward(reg)), c(40L,56L,1L,1L,2L))
    }
})
