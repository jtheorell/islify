context("islify")
data(negImage)
data(posImage)

res <- islify(
    imgDirs = list(negImage, posImage),
    imgNames = c("Neg", "Pos"),
    frameNumFocus = 1,
    sizeCutoff = 40,
    intensityCutoffFocus = 0.035,
    diagnoImgs = FALSE
)

test_that("islify expected output without ref", {
    expect_true(res$fractionOfAll_focus[1] == 0 && 
                    round(res$fractionOfAll_focus[2], 2) == 0.22)
})

res <- islify(
    imgDirs = list(negImage, posImage),
    imgNames = c("Neg", "Pos"),
    frameNumFocus = 1,
    frameNumReference = 2,
    sizeCutoff = 40,
    intensityCutoffFocus = 0.035,
    intensityCutoffReference= 0.05,
    diagnoImgs = FALSE
)

test_that("islify expected output with ref", {
    expect_true(round(res$fractionOfAll_ref[1],2) == 0.22 && 
                    round(res$fractionOfAll_ref[2], 2) == 0.06 && 
                    res$fractionOfAll_focus[1] == 0 && 
                    round(res$fractionOfAll_focus[2], 2) == 0.02 &&
                    res$fractionOfRef_focus[1] == 0 && 
                    round(res$fractionOfRef_focus[2], 2) == 0.39)
})