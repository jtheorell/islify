#' Identification of an intensity cutoff for a whole experiment
#' 
#' With this function, the otherwise integrated intensity cutoff is 
#' externalised. This is only useful in situations where the original values
#' from the imaging have not been perturbed, so that the individual frames have
#' different max values, and are not all compressed to a range from 0 to 1. NB!
#' This requirement is not fulfilled if the files are exported as RGB compliant,
#' as the RGB format always ranges from 0 to 1. 
#' @param imgDirs A vector or list of pathways, including filenames, to the
#' images to be analysed, e.g. "Raw_images/Positive_ctrl.nd2". Formats that are 
#' currently supported are nd2, czi, tiff, png or lists of images, in the form 
#' of three-dimensional arrays, where each layer in the third dimension 
#' represents a color. nd2, czi or non-normalised, integer TIFF are clearly 
#' preferable for memory and resolution purposes. 
#' @param frameNum This identifies which of the frames in the file that 
#' contains the information about the autoantibody binding.
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. 
#' @return An intensity cutoff value.
#' @examples
#' #Load example data and run the function: 
#' data(negImage)
#' data(posImage)
#' getIntensityCutoff(imgDirs = list(negImage, posImage), frameNum = 1)
#' @export getIntensityCutoff
getIntensityCutoff <- function(imgDirs, frameNum, 
                               numOfImgs = "All", numPix = "All"){
  #First, we deal with the special case where the imgDir is not a directory
  #or a list of files, but an individual file, either containing only one
  #color matrix or a list of color matrices. 
  if(is.matrix(imgDirs)){
    imgDirs <- list(imgDirs)
  }
  #In the case where there was only one color matrix in the first place, or
  #where the file is made up of a list of color matrices, we go on to deal with
  #these here. 
  if(is.matrix(imgDirs[[1]])){
    imgDirs <- list(imgDirs)
  }
  median(unlist(lapply(imgDirs, function(x){
    locFile <- importFile(x, frameNum, numOfImgs)
    if(is.numeric(numPix)){
      locFile <- pixReduction(locFile, numPix)
    }
    if(is.list(locFile)){
      locFile <- locFile[[frameNum]]
    }
    auto_thresh(locFile, "tri")
  })))
}