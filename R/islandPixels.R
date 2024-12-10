#Here, we simply count the number of pixels in each of the islands
islandPixels <- function(islandPicture){
  smIsl <- as.data.frame(summary(islandPicture))
  colnames(smIsl) <- c("row", "column", "value")
  splitIslands <- split(smIsl, f = smIsl$value)
  lapply(splitIslands, nrow)
}