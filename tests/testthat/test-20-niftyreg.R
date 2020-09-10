context("Registration")

test_that("Core registration functions work", {
    options(RNiftyReg.threads=2L)

    if (system.file(package="loder") == "")
        skip("The \"loder\" package is not available")
    else
    {
        house <- loder::readPng(system.file("extdata","house.png",package="RNiftyReg"))
        affine <- buildAffine(skews=0.1, source=house, target=house)
        skewedHouse <- applyTransform(affine, house)
        
        reg <- niftyreg(skewedHouse, house, symmetric=FALSE)
        expect_equal(forward(reg)[1,2], 0.1, tolerance=0.1)
        
        reg <- niftyreg(skewedHouse, house, symmetric=TRUE)
        expect_equal(forward(reg)[1,2], 0.1, tolerance=0.1)
        
        # Only true for the symmetric case
        expect_equivalent(invertAffine(forward(reg)), reverse(reg), tolerance=0.0001)
        
        # Hopefully registration has improved the NMI!
        expect_lt(similarity(skewedHouse,house), similarity(reg$image,house))
        
        skip_on_cran()
        skip_on_travis()
        
        reg <- niftyreg(skewedHouse, house, scope="nonlinear", init=forward(reg))
        expect_equal(dim(forward(reg)), c(47L,59L,1L,1L,2L))
    }
})
