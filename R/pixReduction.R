pixReduction <- function(locFile, numPix) {
    rangeList <- lapply(c(nrow, ncol), function(x) {
        topVal <- sample(x(locFile), 1)
        if (topVal < numPix) {
            locVals <- seq(1, numPix)
        } else {
            locVals <- seq((topVal - (numPix - 1)), topVal)
        }
    })
    if (is.matrix(locFile)) {
        locFile[rangeList[[1]], rangeList[[2]]]
    } else {
        locFile[rangeList[[1]], rangeList[[2]], ]
    }
}
