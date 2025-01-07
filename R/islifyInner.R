islifyInner <- function(imgNum, imgDirs, frameNum, sizeCutoff, 
                              diagnoImgs, truncTo, outDir, imgNames, numPix, 
                              highNoise, intensityCutoff, fraction,
                              numOfImgs, reportIntensity, otherPlotDat = FALSE,
                              returnMat, fromIslifyOuter = FALSE){
  imgDir <- imgDirs[[imgNum]]
  locFile <- importFile(imgDir, frameNum, numOfImgs, 
                        fromIslifyOuter = fromIslifyOuter)
  #Now, in cases where the data already has been compressed, we need to 
  #untegerize it for the functions below to work. 
  if(max(locFile) <= 1){
    rescaledList <- rescaleFrom01(locFile, intensityCutoff)
    locFile <- rescaledList[[1]]
    intensityCutoff <- rescaledList[[2]]
  }
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
  if(highNoise){
    intensityCutoff2 <- auto_thresh(round(locFileClean*100), 
                                    "tri", ignore_black = TRUE)/100
    locFileClean[which(locFile <= intensityCutoff2)] <- 0 
  }
  
  locFile01 <- locFileClean
  resList <- list()
  #We also include information about the used filter threshold(s)
  resList$intensityCutoff <- intensityCutoff
  if(highNoise){
    resList$secondIntensityCutoffRed <- intensityCutoff2
  }
  
  #In the case that there is nothing but noise below background, a negative
  #result is thrown back here
  if(length(which(locFile01 != 0)) == 0){
    resList$positivePixels <- 0
    if(reportIntensity){
      resList$intensitySum <- 0
    }
  } else {
    locFile01[which(locFile01 != 0)] <- 1
    locFileIslified <- islandNaming(locFile01)
    
    #Importantly, we also need to remove the now zero points from the 
    #locFileClean, as some points are excluded by the dbscan function. 
    locFileClean[which(locFileIslified == 0)] <- 0
    locFileIslifiedSparse <- Matrix(locFileIslified, sparse = TRUE)
    diameters <- islandDiameter(locFileIslifiedSparse, statistic = "min")
    if(any(diameters > sizeCutoff)){
      locFileBigIsles <- islandRemoval(locFileIslifiedSparse, diameters,
                                       sizeCutoff,
                                       origDims = dim(locFileClean))
      if(fromIslifyOuter == "first"){
        #Now, we calculate the number of pixels per island
        bigIslandPixels <- islandPixels(locFileBigIsles)
        
        ##Here, we are just going to sum the number of pixels that are positive.
        resList$positivePixels <- sum(unlist(bigIslandPixels))
      } else {
        #Here come a few steps that are not relevant if we are running the 
        #outer function for the co-localised protein. 
        
        #Now, we introduce two new statistics, namely the filling 
        #of the islands as well as the removal of all pixels below a 
        #threshold within each island. After this, we will exclude 
        #islands that are fully filled and homogeneous.
        locFileBigIslesFull <- islandFill(locFileBigIsles, 
                                          origDims = dim(locFileClean))
        
        #Now, how many pixels do we get now?
        bigIslandFullPixels <- islandPixels(locFileBigIslesFull)
        
        locFileClean[which(locFileBigIsles == 0)] <- 0
        locHighIsles <- islandPeaks(locFileBigIsles,
                                    Matrix(locFileClean)@x,
                                    fraction)
        diametersHigh <- islandDiameter(locHighIsles, statistic = "min")
        highIslePixels <- islandPixels(locHighIsles)   
        
        #Now, we are going to create an interesting criterion: only islands 
        #that retain their diameter even after all values below the fraction
        #of the top value are kept. 
        
        locFileBigIsles <- islandRemoval(locFileBigIsles, 
                                         diametersHigh,
                                         sizeCutoff,
                                         origDims = dim(locFileClean))
        
        if(reportIntensity){
          #Here, we add information about the sum of 
          #intensity in the big islands
          locFileBigIslesMat <- as.matrix(locFileBigIsles)
          locFileClean[which(locFileBigIslesMat == 0)] <- 0
          resList$intensitySum <- sum(locFileClean)
        }
      }
    } else {
      resList$positivePixels <- 0
      resList$islandDensity <- 1
      if(reportIntensity){
        resList$intensitySum <- 0
      }
    }

      if(resList$positivePixels > 0){
        bigClean <- as.matrix(locFileClean)
        locFileBigIslesMat <- as.matrix(locFileBigIsles)
        bigClean[which(locFileBigIslesMat == 0)] <- 0
      } else {
        bigClean <- matrix(0, nrow(locFileClean), ncol(locFileClean))
      }
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
    
    locFile01 <- plotDat
    if(truncTo == "max"){
      truncVal <- max(plotDat)
    } else {
      truncVal <- truncTo
      plotDat[which(plotDat > truncVal)] <- truncVal
    }
    
    locFile01 <- plotDat/truncVal
    #And here, we invert the colors.
    locFile10 <- abs(locFile01-1)
    plotFile[,,1] <- plotFile[,,2] <- plotFile[,,3] <- locFile10
    
    #And now to the magic: we will remove all values from the green and
    #blue frames, if there are any selected cells. Else, this step is 
    #excluded. 
    if(resList$positivePixels > 0){
      #Now, we start by not flattening, but still decreasing the range of
      #the positive values. 
      locFile0.4 <- (locFile01*0.6)+0.4
      plotFile[,,1][which(bigClean != 0)] <- locFile0.4[which(bigClean != 0)]
      plotFile[,,2][which(bigClean != 0)] <- 0
      plotFile[,,3] <- plotFile[,,2]
    } else {
      plotFile <- locFile0.4
    }
    writePNG(plotFile,paste0(outDir, "/", imgNames[imgNum], ".png"))
  }
  if(returnMat){
    list(unlist(resList), locFileClean)
  } else {
    unlist(resList)
  }
}