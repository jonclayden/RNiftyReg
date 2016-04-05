context("Registration")

test_that("Core registration functions work", {
    if (system.file(package="png") == "")
        skip("The \"png\" package is not available")
    else
    {
        house <- png::readPNG(system.file("extdata","house.png",package="RNiftyReg"))
        affine <- buildAffine(skews=0.1, source=house, target=house)
        skewedHouse <- applyTransform(affine, house)
        
        reg <- niftyreg(skewedHouse, house, symmetric=FALSE)
        expect_that(forward(reg)[1,2], equals(0.1,tolerance=0.1))
        
        reg <- niftyreg(skewedHouse, house, symmetric=TRUE)
        expect_that(forward(reg)[1,2], equals(0.1,tolerance=0.1))
        
        # Only true for the symmetric case
        expect_that(invertAffine(forward(reg)), equals(reverse(reg),tolerance=0.0001,check.attributes=FALSE))
        
        # Hopefully registration has improved the NMI!
        expect_that(similarity(skewedHouse,house), is_less_than(similarity(reg$image,house)))
        
        skip_on_cran()
        
        reg <- niftyreg(skewedHouse, house, scope="nonlinear", init=forward(reg))
        expect_that(dim(forward(reg)), equals(c(47L,59L,1L,1L,2L)))
    }
})
