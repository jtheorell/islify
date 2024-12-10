#Median island intensity is calculated here
islandIntensity <- function(picture, islandPicture){
  if(length(picture@x) != length(islandPicture@x)){
    stop("The picture and he island picture differ in size")
  }
  smRaw <- as.data.frame(summary(picture))
  smIsl <- as.data.frame(summary(islandPicture))
  colnames(smRaw) <- colnames(smIsl) <- c("row", "column", "value")
  splitRaw <- split(smRaw, f = smIsl$value)
  unlist(lapply(split(smRaw$value, f = smIsl$value), median))
}