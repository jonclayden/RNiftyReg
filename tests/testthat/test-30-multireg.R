context("Multiple registration")

test_that("Multiple registration works", {
    if (system.file(package="png") == "")
        skip("The \"png\" package is not available")
    else
    {
        house <- png::readPNG(system.file("extdata","house.png",package="RNiftyReg"))
        affine <- buildAffine(skews=0.1, source=house, target=house)
        skewedHouse <- applyTransform(affine, house)
        
        colourHouse <- png::readPNG(system.file("extdata","house_colour.png",package="RNiftyReg"))
        skewedColourHouse <- applyTransform(affine, colourHouse)
        
        # Rec. 709 luma RGB-to-greyscale coefficients (as used by ImageMagick)
        expect_that(skewedHouse[66,76], equals(weighted.mean(skewedColourHouse[66,76,],c(0.2126,0.7152,0.0722)), tolerance=0.05))
        
        skip_on_os("solaris")
        
        reg <- niftyreg(skewedColourHouse, house, symmetric=FALSE)
        expect_that(forward(reg,2)[1,2], equals(0.1,tolerance=0.05))
        
        skip_on_cran()
        
        reg <- niftyreg(skewedColourHouse, house, scope="nonlinear", init=forward(reg), symmetric=FALSE, maxIterations=300)
        expect_that(dim(forward(reg,2)), equals(c(40L,56L,1L,1L,2L)))
    }
})
