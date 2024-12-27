#' Identification of a suitable size cutoff for cells
#' 
#' This function is needed in situations where the cell size is not 
#' known. What it does is that it identifies the median cell nucleus size
#' in a provided file. In this case, therefore, it is important that the frame
#' number identifies a frame with nuclear staining, such as DAPI or similar. 
#' @param imgDirs A vector or list of pathways, including filenames, to the
#' images to be analysed, e.g. "Raw_images/Positive_ctrl.nd2". Formats that are 
#' currently supported are nd2, czi, tiff, png or lists of images, in the form 
#' of three-dimensional arrays, where each layer in the third dimension 
#' represents a color. nd2, czi or non-normalised, integer TIFF are clearly 
#' preferable for memory and resolution purposes. 
#' @param frameNum This identifies which of the frames in the file that 
#' contains the information about the nuclear staining. In RGB-compliant files, 
#' this is often number 3. 
#' @param noiseThreshold This decides under what size in pixels that an islet is
#' considered noise. Defaults to 10 pixels.
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. 
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @return A size cutoff value in pixels, which corresponds to the median of the
#' island diameters above the noise threshold. 
#' @examples
#' #Load example data and run the function: 
#' data(negImage)
#' getSizeCutoff(imgDirs = list(negImage), frameNum = 3)
#' @export getSizeCutoff
getSizeCutoff <- function(imgDirs, frameNum, noiseThreshold = 10, 
                          numPix = "All", numOfImgs = "All"){
    median(unlist(lapply(imgDirs, sizeCutoffGenInner, frameNum, noiseThreshold, 
                         numPix, numOfImgs)))
}

sizeCutoffGenInner <- function(imgDir, frameNum, noiseThreshold, numPix, 
                               numOfImgs){
  locFile <- importFile(imgDir, frameNum, numOfImgs)
  #Now, in cases where the data already has been compressed, we need to 
  #untegerize it for the functions below to work. 
  if(max(locFile) <= 1){
    rescaledList <- rescaleFrom01(locFile, intensityCutoff = FALSE)
    locFile <- rescaledList[[1]]
  }
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
    #Here, we remove the vast majority of islands, that are tiny.
    realDiameters <- diameters[-which(diameters <= noiseThreshold)]
    #And of the remaining ones, we pick the median. 
    median(realDiameters)
}
