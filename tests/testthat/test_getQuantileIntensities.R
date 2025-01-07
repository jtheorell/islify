context("getQuantileIntensities")
data(negImage)
res <- round(getQuantileIntensities(list(negImage))[[1]], 2)[,1]
rightAnswer <- c(0.02, 0.47, 0.32)
test_that("getQuantileIntensities expected output", {
  expect_true(all(res == rightAnswer))
})
