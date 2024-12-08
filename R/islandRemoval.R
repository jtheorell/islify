#Here is a function to remove small isles
islandRemoval <- function(islandPicture, islandDiameters, sizeCutoff, origDims){
  smallIslands <- islandDiameters[which(islandDiameters <= sizeCutoff)]
  sm <- as.data.frame(summary(islandPicture))
  colnames(sm) <- c("row", "column", "value")
  
  sm$value[which(sm$value %in% as.numeric(names(smallIslands)))] <- 0
  #Here, we remove all new zeros. 
  smSmall <- sm[which(sm$value > 0),]
  locSM <- sparseMatrix(i = smSmall$row, 
                        j = smSmall$column, 
                        x = smSmall$value,
                        dims = origDims)
}
