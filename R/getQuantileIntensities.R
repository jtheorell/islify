#' Get information about the quantiles for the intensity values for all colors
#' in an image
#' 
#' With this function, the user can get an overview of how the values for the 
#' intensities in a file or set of files are distributed. This is useful mainly
#' to set reasonable cutoff values for the saveImage function, but might also be
#' good for sanity checking the getIntensityCutoff results. 
#' @param imgDirs A vector or list of pathways, including filenames, to the
#' images to be analysed, e.g. "Raw_images/Positive_ctrl.nd2". Formats that are 
#' currently supported are nd2, czi, tiff, png or lists of images, in the form 
#' of three-dimensional arrays, where each layer in the third dimension 
#' represents a color. nd2, czi or non-normalised, integer TIFF are clearly 
#' preferable for memory and resolution purposes. 
#' @param quantiles Which quantiles should be returned? Values between and 
#' including 0 and 1 are accepted. 
#' @return A list of matrices with one row per color in the data and with 
#' thirteen columns of intensities from percentile 0 through to percentile 100. 
#' @examples
#' #Load example data and run the function: 
#' data(negImage)
#' getQuantileIntensities(list(negImage))
#' 
#' @export getQuantileIntensities
getQuantileIntensities <- function(imgDirs, quantiles = c(0.99)){
  lapply(imgDirs, function(x){
    locFile <- importFile(x, frameNum = "All", numOfImgs = "All")
    locMat <- do.call("rbind", lapply(seq(1, dim(locFile)[3]), function(y){
      quantile(x[,,y], quantiles)
    }))
    locDf <- as.data.frame(locMat)
    colnames(locDf) <- paste0("Percent_", gsub("|%", "", colnames(locMat)))
    locDf
  })
}
