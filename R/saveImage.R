#' Saving individually coloured png images from raw image input
#' 
#' THis function simply generates pngs from raw image files, such as .nd2.
#' 
#' @param imgDirs The directory containing the image files. NB! Only nd2
#' files are currently supported.
#' @param imgNames The names of the images. This argument is used both for
#' naming of the rows in the output matrix, but also for naming of the images
#' if such are generated using the diagnoImgs flag. 
#' @param frameNums This identifies which of the frames in the files that 
#' should be plotted. Defaults to "All", which then also generates a merge. 
#' @param frameCols Which color does the frame/frames in question have? 
#' Strangely, in some file formats such as nd2, it varies - not all are 
#' "R-G-B", but "G-R-B" also occurs. Can take the values "R"/"Red", "G"/"Green",
#' and "B"/"Blue". Should be a vector the same length as frameNums. If nothing
#' is provided, it defaults to only saving the merge. 
#' @param normTo This argument has three possible values: "bit", "max" or a
#' numeric value/a numeric vector the length of frameNums. Max is png standard,
#' whereas "bit" makes all pictures in a set comparable, but will make it hard
#' to see the subtle shifts often present in negative files. The third option
#' can be used if the max value in a set of files should be used, but this then
#' needs to be externally defined. If this value is set to a lower value than
#' the true max, the data will be truncated and a warning thrown. 
#' @param numBit How many bit are the nd2 files? Can be found using 
#' \code{\link{RBioFormats::read.metadata}} function, in the 
#' "$coreMetadata$bitsPerPixel" slot. NB! If more than one image is present in
#' the nd2 file, one needs to refer to the first slot in the list, i.e., sub-
#' setting the metadata objects with double hard brackets.
#' @param outDir The directory that the images should be saved in. Subfolders
#' are created within this directory for each color. 
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. NB! If this command is used, then the images
#' will not be identical from round to round, and will not be comparable to the
#' output of other functions, as a random subset of the pictures are used. 
#' @return Named png images that are scaled to the max allowed value in the 
#' raw image file and coloured according to frameCol.
#' 
#' @export saveImage
saveImage <- function(imgDirs, imgNames, frameNums = "All", frameCols,
                      normTo = "bit", 
                      numBit = 16, outDir, numOfImgs = "All", 
                      numPix = "All"){
  if(frameNums == "All"){
    for(i in c("Merged_colors", frameCols)){
      dir.create(paste0(outDir, "/",i))
    }
    
  }
  lapply(1:length(imgDirs), function(x){
    imgDir <- imgDirs[[x]]
    imgName <- imgNames[[x]]
    #We start by the default, i.e. that we want to plot everything. 
    if(frameNums == "All"){
      locFile <- importFile(imgDir, 
                            frameNum = "All", 
                            numOfImgs = numOfImgs)
      #Here, if max is 1, all other normalisation information will be discarded,
      #with a warning. 
      if(max(locFile) == 1 && (normTo[1] == "bit" || normTo[1] > 1)){
        message("This file has a max value of 1, and therefore this will be",
        "used for truncation purposes, if the provided ")
        normTo <- 1
      }
      
      if(is.numeric(numPix)){
        locFile <- pixReduction(locFile, numPix)
      }
      if(length(normTo) == 1){
        if(normTo == "bit"){
          normTo <- rep(2^numBit, length(frameCols))
        } else if (normTo == "max"){
          normTo <- vapply(seq_along(frameCols), function(y) max(locFile[,,y]),1)
        } else{
          normTo <- rep(normTo, length(frameCols))
        }
      }
      if(length(normTo) != length(frameCols)){
        stop("Mismatch between the length of frameCols and normTo. Change this.")
      }
      #First, we make the plots for individual colors
      for(i in seq_along(frameCols)){
        frameCol <- frameCols[i]
        saveImageInner(locImg = locFile[,,i], 
                       frameCol = frameCol, 
                       normToInner = normTo[i], 
                       imgName = imgName, 
                       outDir = outDir)
      }
      #Then the common plot. 
      #Here, we normalise to the range of the possible max value.
      locImgNorm <- locFile
      for(i in seq_along(frameCols)){
        locImgNorm[,,i] <- locFile[,,i]/normTo[i]
      }
      
      frameNums <- vapply(frameCols, colToNum, 1)
      
      #Now, if a normTo value lower than the max is used, a warning is thrown and
      #the data is truncated
      if(max(locImgNorm) > 1){
        warning("The max value in the file exceeds the normTo value and is truncated")
        locImgNorm[which(locImgNorm > 1)] <- 1
      }
      #And now, we create a full rgb array. 
      locArray <- array(0,dim= c(nrow(locImgNorm),
                                 ncol(locImgNorm),
                                 3))
      for(i in seq_along(frameCols)){
        locArray[,,frameNums[i]] <- locImgNorm[,,i]
      }

      writePNG(locArray,paste0(outDir, "/Merged_colors/", imgName, "_merge.png"))
    } else {
      locFile <- importFile(imgDir, 
                            frameNum = frameNums, 
                            numOfImgs = numOfImgs)
      saveImageInner(locImg = locFile, 
                     frameCol = frameCol, 
                     normToInner = normTo, 
                     imgName = imgName, 
                     outDir = outDir)
    }
    
    
  })
}

saveImageInner <- function(locImg, frameCol, normToInner, imgName, outDir){
  
  if(normToInner == "bit"){
    normToInner <- 2^numBit
  } else if (normToInner == "max"){
    normToInner <- max(locImg)
  }
  #Here, we normalise to the range of the possible max value.
  locImgNorm <- locImg/normToInner
  
  frameNum <- colToNum(frameCol)
  
  #Now, if a normTo value lower than the max is used, a warning is thrown and
  #the data is truncated
  if(max(locImgNorm) > 1){
    warning("The max value in the file exceeds the normTo value and is truncated")
    locImgNorm[which(locImgNorm > 1)] <- 1
  }
  #And now, we create a full rgb array. 
  locArray <- array(0,dim= c(nrow(locImgNorm),
                             ncol(locImgNorm),
                             3))
  
  locArray[,,frameNum] <- locImgNorm
  writePNG(locArray,paste0(outDir, "/",frameCol, 
                           "/", imgName, "_", 
                           frameCol, ".png"))
}

colToNum <- function(frameCol){
  switch(frameCol, 
         "R" = 1, 
         "Red" = 1,
         "red" = 1,
         "G" = 2,
         "Green" = 2,
         "green" = 2,
         "B" = 3,
         "Blue" = 3,
         "blue" = 3)
}