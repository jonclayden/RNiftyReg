# s <- readImage("source.nii.gz")
# t <- readImage("target.nii.gz")
# a <- readImage("alternate.nii.gz")
#
# r1 <- niftyreg(s, t, scope="affine")
# r2 <- niftyreg(s, t, scope="nonlinear", init=forward(r1))
#
# ar1 <- applyTransform(reverse(r1), a)
# ar2 <- applyTransform(reverse(r2), a)
#
# p <- applyTransform(forward(r1), c(32,45,12))
#
# aff <- readAffine("affine.txt")
# ar3 <- applyTransform(aff, a)
