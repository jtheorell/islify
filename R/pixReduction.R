pixReduction <- function(locFile, numPix){
  rangeList <- lapply(c(nrow, ncol), function(x){
    topVal <- sample(x(locFile), 1)
    if(topVal < numPix){
      locVals <- 1:numPix
    } else {
      locVals <- (topVal-(numPix-1)):topVal
    }
  })
  locFile[rangeList[[1]],rangeList[[2]]]
}