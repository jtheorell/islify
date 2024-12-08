importFile <- function(imgDir, frameNum, numOfImgs){
  if(is.list(imgDir)){
    message("This object is interpreted as an image file and not a directory")
    locFile <- imgDir[[frameNum]]
  } else if(grepl("\\.nd2", imgDir)){
    message(imgDir)
    #For this, we need to be clever, as many files can be saved in the same
    #object. For now, we will allow all to go through. 
    locFileSet <- read.image(imgDir, normalize = FALSE)
    if(numOfImgs == "All"){
      numOfImgs <- length(locFileSet@.Data)
    }
    locFile <- do.call("cbind", lapply(seq(1, numOfImgs), function(x){
      locFileSet[[x]]@.Data[,,frameNum]
    }))
  } else if(grepl(".png", imgDir)){
    message(imgDir)
    locFile <- readPNG(imgDir)[,,frameNum]
  } else if(grepl(".tif", imgDir)){
    message(imgDir)
    locFile <- readTIFF(imgDir, as.is = TRUE, all = TRUE)[[frameNum]]
    #Here, we handle the situation where the imported file follow RGB standards
    #and thus contains three colors for each 
    if(is.array(locFile)){
      message(paste0("NB! This file is interpreted as an RGB-compliant",
                     " 3-dimensional array. If the number of colors in the",
                     " experiment is not three, the frame number used here",
                     " might mislead the analysis. Be skeptical to the outome", 
                     " in such cases."))
      locFile <- locFile[,,frameNum]
    }
  }
}