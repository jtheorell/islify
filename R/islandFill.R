
#With this function, we are going column by column for each cluster and 
#filling all pixels that lie between the edges. 
islandFill <- function(islandPicture, origDims){
  smIsl <- as.data.frame(summary(islandPicture))
  colnames(smIsl) <- c("row", "column", "island")
  splitIslands <- split(smIsl, f = smIsl$island)
  fullPicture <- do.call("rbind", lapply(splitIslands, function(x){
    locIsland <- x$island[1]
    locSplit <- split(x, x$column)
    fullIsland <- do.call("rbind", lapply(locSplit, function(y){
      minVal <- min(y$row)
      maxVal <- max(y$row)
      locSeq <- seq(minVal, maxVal)
      data.frame("row" = locSeq, "column" = y$column[1])
    }))
    fullIsland$island <- locIsland
    fullIsland
  }))
  
  locSM <- sparseMatrix(i = fullPicture$row, 
                        j = fullPicture$column, 
                        x = fullPicture$island,
                        dims = origDims)
}