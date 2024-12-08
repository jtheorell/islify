#' Saving individually coloured png images from raw image input
#' 
#' THis function simply generates pngs from raw image files, such as .nd2.
#' 
#' @param imgDirs The directory containing the image files. NB! Only nd2
#' files are currently supported.
#' @param imgNames The names of the images. This argument is used both for
#' naming of the rows in the output matrix, but also for naming of the images
#' if such are generated using the diagnoImgs flag. 
#' @param frameNum This identifies which of the frames in the files that 
#' contains the information about the autoantibody binding.
#' @param frameCol Which color does the frame in question have? Strangely, it
#' varies in the nd2 files - not all are "R-G-B", but "G-R-B" also occurs. 
#' Can take the values "R", "G", and "B" for red, green and blue. 
#' @param normTo This argument has three possible values: "bit", "max" or a
#' numeric value. Max is png standard, whereas the former makes all 
#' pictures in a set comparable, but will make it hard to see the subtle shifts
#' often present in negative files. The third option can be used if the max
#' value in a set of files should be used, but this then needs to be externally
#' defined. If this value is set to a lower value than the true max, the data
#' will be truncated and a warning thrown. 
#' @param numBit How many bit are the nd2 files? Can be found using 
#' \code{\link{RBioFormats::read.metadata}} function, in the 
#' "$coreMetadata$bitsPerPixel" slot. NB! If more than one image is present in
#' the nd2 file, one needs to refer to the first slot in the list, i.e., sub-
#' setting the metadata objeects with double hard brackets.
#' @param outDir The directory that the diagnoImages should be saved in. Only
#' used together with diagnoImgs = TRUE.
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. NB! If this command is used, then the images
#' will not be identical from round to round, and will not be comparable to the
#' output of other functions, as a randdom subset of the pictures are used. 
#' @return Named png images that are scaled to the max allowed value in the 
#' raw image file and coloured according to frameNum (1=red, 2=green, 3=blue).
#' @export saveImage
saveImage <- function(imgDirs, imgNames, frameNum, frameCol,
                      normTo = "bit", 
                      numBit = 16, outDir, numOfImgs = "All", 
                      numPix = "All"){
  lapply(1:length(imgDirs), function(x){
    imgDir <- imgDirs[[x]]
    locFile <- importFile(imgDir, frameNum, numOfImgs)
    if(is.numeric(numPix)){
      locFile <- pixReduction(locFile, numPix)
    }
    if(normTo == "bit"){
      normTo <- 2^numBit
    } else if (normTo == "max"){
      normTo <- max(locFile)
    }
    
    #Here, we normalise to the range of the possible max value.
    locFileNorm <- locFile/normTo
    
    #Now, if a normTo value lower than the max is used, a warning is thrown and
    #the data is truncated
    if(max(locFileNorm) > 1){
      warning("The max value in the file exceeds the normTo value and is truncated")
      locFileNorm[which(locFileNorm > 1)] <- 1
    }
    #And now, we create a full rgb array. 
    locArray <- array(0,dim= c(nrow(locFileNorm),
                               ncol(locFileNorm),
                               3))
    colNum <- switch(frameCol, 
                     "R" = 1, 
                     "G" = 2,
                     "B" = 3)
    locArray[,,colNum] <- locFileNorm
    
    writePNG(locArray,paste0(outDir, "/", imgNames[x], ".png"))
  })
}