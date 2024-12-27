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
#' @return A matrix with one row per color in the data and with thirteen columns
#' of intensities from quantile 0 through to quantile 100. 
#' @examples
#' #Load example data and run the function: 
#' data(negImage)
#' getImageIntensities(list(negImage))
#' 
#' @export getImageIntensities
getImageIntensities <- function(imgDirs){
  lapply(imgDirs, function(x){
    locFile <- importFile(x, frameNum = "All", numOfImgs = "All")
    do.call("rbind", lapply(seq(1, dim(locFile)[3]), function(y){
      quantile(x[,,y], c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1))
    }))
  })
}
