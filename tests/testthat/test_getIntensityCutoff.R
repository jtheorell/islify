context("getIntensityCutoff")
data(negImage)
res <- getIntensityCutoff(imgDirs = list(negImage), frameNum = 2)

test_that("getIntensityCutoff expected output", {
    expect_true(res == 0.05)
})

res <- getIntensityCutoff(imgDirs = list(negImage), frameNum = 2,
                          ignore_white = TRUE)

test_that("getIntensityCutoff expected output ignore white TRUE", {
    expect_true(res == 0.05)
})

res <- getIntensityCutoff(imgDirs = list(negImage), frameNum = 2,
                          ignore_white = 0.1)

test_that("getIntensityCutoff expected output ignore white 0.1", {
    expect_true(res == 0.024)
})