islandNaming <- function(picture){
  indices <- which(picture == 1, arr.ind = TRUE)
  locIsles <- dbscan(indices, eps = 1, minPts = 5)
  pictureResult <- picture
  pictureResult[which(picture == 1)] <- locIsles$cluster
  pictureResult
}