context("getSizeCutoff")
data(negImage)
res <- getSizeCutoff(imgDirs = list(negImage), frameNum = 3)

test_that("getSizeCutoff expected output", {
    expect_true(res == 38.5)
})

