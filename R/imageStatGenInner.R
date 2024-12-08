imageStatGenInner <- function(imgNum, imgDirs, frameNum, sizeCutoff, 
                              diagnoImgs, outDir, imgNames, numPix, fixed, 
                              intensityCutoff, numOfImgs, otherPlotDat = FALSE,
                              returnMat){
  imgDir <- imgDirs[[imgNum]]
  locFile <- importFile(imgDir, frameNum, numOfImgs)
  #Here, we restrict the number of pixels, in applicable cases, to a randomly
  #located subset of the frame
  if(is.numeric(numPix)){
    locFile <- pixReduction(locFile, numPix)
  }
  locFileClean <- locFile
  
  if(is.logical(intensityCutoff) && intensityCutoff == TRUE){
    #This specific triangle method seems to work the best
    intensityCutoff <- auto_thresh(locFile, "tri")
  }
  
  locFileClean[which(locFile <= intensityCutoff)] <- 0
  #Now, we are going to introduce an optional second filtering step, where 
  #we disregard the information from the areas devoid of cells, and also remove
  #the noise from the non-important cells. 
  if(fixed){
    intensityCutoff2 <- auto_thresh(round(locFileClean*100), 
                                    "tri", ignore_black = TRUE)/100
    locFileClean[which(locFile <= intensityCutoff2)] <- 0 
  }
  
  locFile01 <- locFileClean
  
  #In the case that there is nothing but noise below background, a negative
  #result is thrown back here
  if(length(which(locFile01 != 0)) == 0){
    negOutput()
  } else {
    locFile01[which(locFile01 != 0)] <- 1
    locFileIslified <- islandNaming(locFile01)
    
    #Importantly, we also need to remove the now zero points from the 
    #locFileClean, as some points are exclded by the dbscan function. 
    locFileClean[which(locFileIslified == 0)] <- 0
    locFileIslifiedSparse <- Matrix(locFileIslified, sparse = TRUE)
    diameters <- islandDiameter(locFileIslifiedSparse)
    if(any(diameters > sizeCutoff)){
      locFileBigIsles <- islandRemoval(locFileIslifiedSparse, diameters,
                                       sizeCutoff,
                                       origDims = dim(locFileClean))
      bigIslandIndexes <- unique(locFileBigIsles@x)
    } else {
      locFileBigIsles <- new("nullDgCMatrix")
      bigIslandIndexes <- 0
    }

      resList <- list()
      resList$nIsles <- length(unique(locFileBigIsles@x))
      
      #resList$intensitySum <- sum(locFileClean)
      #
     
      if(resList$nIsles > 0){
        bigClean <- as.matrix(locFileClean)
        locFileBigIslesMat <- as.matrix(locFileBigIsles)
        bigClean[which(locFileBigIslesMat == 0)] <- 0
      } else {
        bigClean <- matrix(0, nrow(locFileClean), ncol(locFileClean))
      }
      #At this stage, we can, if it is prompted, throw back a picture of the 
      #selected islands in the context of the full file. 
      if(diagnoImgs){
        plotFile <- array(matrix(0,nrow(locFile),ncol(locFile)), 
                          c(nrow(locFile), ncol(locFile), 3))
        #In the case where another file has been added that should be the background
        #plot, i.e. when we are using the imageStatGenOuter, we add this here
        if(is.logical(otherPlotDat) == FALSE){
          plotDat <- otherPlotDat
        } else {
          plotDat <- locFile
        }
        locFile01 <- plotDat/max(plotDat)
        #And here, we invert the colors.
        locFile10 <- abs(locFile01-1)
        plotFile[,,1] <- plotFile[,,2] <- plotFile[,,3] <- locFile10
        #And now to the magic: we will remove all values from the green and
        #blue frames, if there are any selected cells. Else, this step is 
        #excluded. 
        if(resList$nIsles > 0){
          #Now, we start by not flattening, but still decreasing the range of
          #the positive values. 
          locFile0.4 <- (locFile01*0.6)+0.4
          plotFile[,,1][which(bigClean != 0)] <- locFile0.4[which(bigClean != 0)]
          plotFile[,,2][which(bigClean != 0)] <- 0
          plotFile[,,3] <- plotFile[,,2]
        }
        writePNG(plotFile,paste0(outDir, "/", imgNames[imgNum], ".png"))
      }
      if(returnMat){
        list(unlist(resList), locFileClean)
      } else {
        unlist(resList)
      }
      
  }
}