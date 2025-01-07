context("getIntensityCutoff")
data(negImage)
res <- getIntensityCutoff(imgDirs = list(negImage), frameNum = 2)

test_that("getIntensityCutoff expected output", {
    expect_true(res == 0.05)
})

