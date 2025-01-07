#The point of this function is to provide a per-island height filter, as the
#intensities might vary considerably for each cell.
islandPeaks <- function(islandPicture, intensityVec, fraction){
  origDims <- dim(islandPicture)
  smIsl <- as.data.frame(summary(islandPicture))
  colnames(smIsl) <- c("row", "column", "island")
  splitIslands <- split(smIsl, f = smIsl$island)
  splitIntensities <- split(intensityVec, f = smIsl$island)

  redPicture <- 
    do.call("rbind",
            lapply(seq_along(splitIslands), function(x){
              locIsland <- splitIslands[[x]]
              locIntensities <- splitIntensities[[x]]
              maxVal <- max(locIntensities)
              fracVal <- maxVal*fraction
              locIsland[-which(locIntensities < fracVal),]
    }))
locSM <- sparseMatrix(i = redPicture$row, 
                      j = redPicture$column, 
                      x = redPicture$island,
                      dims = origDims)
}