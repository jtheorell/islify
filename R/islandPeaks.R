# The point of this function is to provide a per-island height filter, as the
# intensities might vary considerably for each cell. We use the 99th
# percentile of the cells here, to avoid the very extreme events having an
# unduly influence on the results.
islandPeaks <- function(islandPicture, intensityVec, ringFrac, flatFrac) {
    origDims <- dim(islandPicture)
    smIsl <- as.data.frame(summary(islandPicture))
    colnames(smIsl) <- c("row", "column", "island")
    splitIslands <- split(smIsl, f = smIsl$island)
    splitIntensities <- split(intensityVec, f = smIsl$island)

    redPicture <-
        do.call(
            "rbind",
            lapply(seq_along(splitIslands), function(x) {
                locIsland <- splitIslands[[x]]
                locIntensities <- splitIntensities[[x]]
                maxVal <- quantile(locIntensities, 0.99)
                #Here, we are actually evaluating two different aspects: 
                #First, we check whether the median absolute deviation
                #of the values is very low compared to the range from max
                #to 0. If so, the island is deemed flat and reduced to a 
                #single pixel. 
                locMad <- mad(locIntensities)
                if(locMad < flatFrac*maxVal){
                    locIsland[1,]
                } else {
                    minVal <- min(locIntensities)
                    range <- maxVal - minVal
                    fracVal <- (ringFrac * range) + minVal
                    if(any(locIntensities < fracVal)){
                        locIsland[-which(locIntensities < fracVal), ]
                    } else {
                        locIsland
                    }
                    
                }
            })
        )
    locSM <- sparseMatrix(
        i = redPicture$row,
        j = redPicture$column,
        x = redPicture$island,
        dims = origDims
    )
}
