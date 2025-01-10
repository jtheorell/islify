saveImageInner <- function(locImg, frameCol, truncToInner,
                           imgName, outDir) {
    if (truncToInner == "max") {
        truncToInner <- max(locImg)
    }
    # Here, we normalise to the range of the possible max value.
    locImgNorm <- locImg / truncToInner

    frameNum <- colToNum(frameCol)

    if (max(locImgNorm) > 1) {
        locImgNorm[which(locImgNorm > 1)] <- 1
    }
    # And now, we create a full rgb array.
    locArray <- array(0, dim = c(
        nrow(locImgNorm),
        ncol(locImgNorm),
        3
    ))

    locArray[, , frameNum] <- locImgNorm
    writePNG(locArray, paste0(
        outDir, "/", frameCol,
        "/", imgName, "_",
        frameCol, ".png"
    ))
}
