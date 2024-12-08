#' Identification of an intensity cutoff for a whole experiment
#' 
#' With this function, the otherwise integrated intensity cutoff is 
#' externalised. This is ONLY useful in situations where the original values
#' from the imaging have not been perturbed, so that the individual frames have
#' different max values, and are not all compressed to a range from 0 to 1. NB!
#' This requirement is not fulfilled if the files are exported as RGB compliant,
#' as the RGB format always ranges from 0 to 1. 
#' @param imgDirs The directory where the file(s) that will be used to generate
#' the threshold is kept. If more than one file, then a median value is
#' generated. Currently, this function is only implemented to work with 
#' raw-value non-RGB TIFF files. 
#' There is an option here, which is mainly for example purposes, using the
#' provided dataset, and that is that a list of image files can be provided 
#' directly. These should then be in the form of a list of lists, each sublist 
#' being composed of matrices representing the colors. See example for useage. 
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
#' @export intensityCutoffGen
#' @examples
#' See example in \code{\link{imageStatGen}}
intensityCutoffGen <- function(imgDirs, frameNum, 
                               numOfImgs = "All", numPix = "All"){
  #First, we deal with the special case where the imgDir is not a directory
  #or a list of files, but an individual file, either containing only one
  #color matrix or a list of color matrices. 
  if(is.matrix(imgDirs)){
    imgDirs <- list(imgDirs)
  }
  #In the case where there was only one color matrix in the first place, or
  #where the file is made up of alist of color matrices, we go on to deal with
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