getSizeCutoffInner <- function(imgDir, frameNum, noiseThreshold, numPix,
                               numOfImgs) {
    locFile <- importFile(imgDir, frameNum, numOfImgs)
    # Now, in cases where the data already has been compressed, we need to
    # untegerize it for the functions below to work.
    if (max(locFile) <= 1) {
        rescaledList <- rescaleFrom01(locFile,
            intensityCutoff = FALSE,
            printMessage = FALSE
        )
        locFile <- rescaledList[[1]]
    }
    # Here, we restrict the number of pixels, in applicable cases, to a
    # randomly located subset of the frame
    if (is.numeric(numPix)) {
        locFile <- pixReduction(locFile, numPix)
    }
    intensityCutoff <- auto_thresh(locFile, "tri")
    locFileClean <- locFile
    locFileClean[which(locFile <= intensityCutoff)] <- 0
    locFile01 <- locFileClean
    locFile01[which(locFile01 != 0)] <- 1
    locFileIslified <- Matrix(islandNaming(locFile01), sparse = TRUE)
    diameters <- islandDiameter(locFileIslified)
    # Here, we remove the vast majority of islands, that are tiny.
    realDiameters <- diameters[-which(diameters <= noiseThreshold)]
    # And of the remaining ones, we pick the median.
    median(realDiameters)
}
