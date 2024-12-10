#' Identification of a suitable size cutoff for cells
#' 
#' This function is needed in situations where the cell size is not 
#' known. What it does is that it identifies the median cell nucleus size
#' in a provided file. In this case, therefore, it is important that the frame
#' number identifies a frame with nuclear staining, such as DAPI or similar. 
#' @param imgDirs The directory containing the image file/files that should be
#' used to identify the cell size. In this case, one or a few such files should
#' suffice, as the cells are expected to behave very similarly in all frames. 
#' There is an option here, which is mainly for example purposes, using the
#' provided dataset, and that is that a list of image files can be provided 
#' directly. These should then be in the form of a list of lists, each sublist 
#' being composed of matrices representing the colors. See example for useage. 
#' @param frameNum This identifies which of the frames in the file that 
#' contains the information about the nuclear staining. In RGB files, this is 
#' often number 3. 
#' @param returAllVals Should all the island diameters, including the very tiny,
#' noisy ones be returned, or should a focused value be returned? 
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. 
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @return A size cutoff value in pixels, which corresponds to either:
#'  - the media of the "real" island diameters, i.e. after exclusion 
#'  of noisy, minimal islands of 1-10 pixel diameters.
#'  - a vector of all island diameters, including the noise, if returAllVals
#'  is set to TRUE. 
#' @export sizeCutoffGen
#' @examples
#' See example in \code{\link{imageStatGen}}
sizeCutoffGen <- function(imgDirs, frameNum, returnAllVals = FALSE, 
                          numPix = "All", numOfImgs = "All"){
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
    median(unlist(lapply(imgDirs, sizeCutoffGenInner, frameNum, returnAllVals, 
                         numPix, numOfImgs)))
}

sizeCutoffGenInner <- function(imgDir, frameNum, returnAllVals, numPix, 
                               numOfImgs){
  locFile <- importFile(imgDir, frameNum, numOfImgs)
  #Here, we restrict the number of pixels, in applicable cases, to a randomly
  #located subset of the frame
  if(is.numeric(numPix)){
    locFile <- pixReduction(locFile, numPix)
  }
    #This specific triangle method with the following settings seem to work best
    intensityCutoff <- auto_thresh(locFile, "tri")
    locFileClean <- locFile
    locFileClean[which(locFile <= intensityCutoff)] <- 0
    locFile01 <- locFileClean
    locFile01[which(locFile01 != 0)] <- 1
    locFileIslified <- Matrix(islandNaming(locFile01), sparse = TRUE)
    diameters <- islandDiameter(locFileIslified)
    if(returnAllVals){
      diameters
    } else {
      #Here, we remove the vast majority of islands, that are tiny, 1-10
      #pixel things
      realDiameters <- diameters[-which(diameters <=10)]
      #And of the remaining ones, we pick the median. 
      median(realDiameters)
    }
}
