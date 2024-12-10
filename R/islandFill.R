
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
  
  #This code is more elegant, but alas, slower, so it is commented out for now.
  #This might be converted to c++ code though. 
  #splitIslandCols <- split(smIsl, f = paste(smIsl$column, smIsl$island))
  #fullPicture <- do.call("rbind", lapply(splitIslandCols, function(x){
  #    data.frame("row" = seq(min(x$row), max(x$row)), 
  #               "column" = x$column[1], 
  #               "island" = x$island[1])
  #}))

  locSM <- sparseMatrix(i = fullPicture$row, 
                        j = fullPicture$column, 
                        x = fullPicture$island,
                        dims = origDims)
}