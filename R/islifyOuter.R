islifyOuter <- function(imgNum, imgDirs,
                        frameNumFocus, frameNumReference,
                        sizeCutoff,
                        diagnoImgs, truncTo, outDir,
                        imgNames, numPix, highNoise,
                        intensityCutoffFocus,
                        intensityCutoffReference,
                        threshold_method,
                        ignore_white,
                        fraction,
                        reportIntensity,
                        numOfImgs) {
    # We start by analysing the reference picture
    refImageDat <- islifyInner(
        imgNum = imgNum,
        imgDirs = imgDirs,
        frameNum = frameNumReference,
        sizeCutoff = sizeCutoff,
        diagnoImgs = FALSE,
        truncTo = FALSE,
        numPix = numPix,
        highNoise = highNoise,
        intensityCutoff = intensityCutoffReference,
        threshold_method = threshold_method,
        ignore_white = ignore_white,
        fraction = 1,
        reportIntensity = FALSE,
        fromIslifyOuter = "first",
        numOfImgs = numOfImgs
    )
    # Now, we locally import the focus image
    focImageRaw <- importFile(imgDirs[[imgNum]], frameNumFocus, numOfImgs)
    focImageFiltered <- focImageRaw
    focImageFiltered[which(refImageDat[[2]] == 0)] <- 0
    # Here, we construct a new image with this data.
    focImgFull <- list(abind(focImageFiltered,
        matrix(
            0,
            nrow(focImageFiltered),
            ncol(focImageFiltered)
        ),
        along = 3
    ))
    focImageDat <- islifyInner(
        imgNum = 1,
        imgDirs = focImgFull,
        frameNum = 1,
        sizeCutoff = sizeCutoff,
        diagnoImgs = diagnoImgs,
        truncTo = truncTo,
        outDir = outDir,
        imgNames = imgNames[[imgNum]],
        numPix = numPix,
        highNoise = highNoise,
        intensityCutoff = intensityCutoffFocus,
        threshold_method = threshold_method,
        ignore_white = ignore_white,
        fraction = fraction,
        reportIntensity = reportIntensity,
        fromIslifyOuter = "second",
        numOfImgs = numOfImgs,
        otherPlotDat = focImageRaw
    )
    # Now, we re-formulate the output
    refImageResRaw <- as.numeric(refImageDat[[1]])
    focImageResRaw <- as.numeric(focImageDat)
    resultVec <- c(
        "intensityCutoff_ref" = refImageResRaw[1],
        "intensityCutoff_focus" = focImageResRaw[1],
        "fractionOfAll_ref" = refImageResRaw[2],
        "fractionOfAll_focus" = focImageResRaw[2],
        "fractionOfRef_focus" = focImageResRaw[2] / refImageResRaw[2]
    )
    if (highNoise) {
        resultVec <- c(resultVec,
            "secondIntensityCutoff_ref" = refImageResRaw[3],
            "secondIntensityCutoff_focus" = focImageResRaw[3]
        )
    }
    resultVec
}
