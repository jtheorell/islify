imageStatGenOuter <- function(imgNum, imgDirs, frameNumFocus, frameNumReference,
                              sizeCutoff, 
                              diagnoImgs, outDir, imgNames, numPix, highNoise, 
                              intensityCutoffFocus, 
                              intensityCutoffReference, nuclearToCellQuotient,
                              reportIntensity,
                              numOfImgs){
  #We start by analysing the reference picture
  secondImageDat <- imageStatGenInner(imgNum = imgNum,
                                     imgDirs = imgDirs, 
                                     frameNum = frameNumReference,
                                     sizeCutoff = sizeCutoff, 
                                     diagnoImgs = FALSE, 
                                     numPix = numPix, 
                                     highNoise = highNoise, 
                                     intensityCutoff = intensityCutoffReference,
                                     nuclearToCellQuotient = 
                                       nuclearToCellQuotient,
                                     reportIntensity = FALSE,
                                     numOfImgs = numOfImgs,
                                     returnMat = TRUE)
  #Now, we locally import the focus image
  primaryImageRaw <- importFile(imgDirs[[imgNum]], frameNumFocus, numOfImgs)
  primaryImageFiltered <- primaryImageRaw
  primaryImageFiltered[which(secondImageDat[[2]] == 0)] <- 0
  primaryImageDat <- imageStatGenInner(imgNum = 1,
                                   imgDirs = list(list(primaryImageFiltered)), 
                                   frameNum = 1,
                                   sizeCutoff = sizeCutoff, 
                                   diagnoImgs = diagnoImgs, 
                                   outDir = outDir, 
                                   imgNames = imgNames[[imgNum]],
                                   numPix = numPix, 
                                   highNoise = highNoise, 
                                   intensityCutoff = intensityCutoffFocus,
                                   nuclearToCellQuotient = 
                                     nuclearToCellQuotient,
                                   reportIntensity = reportIntensity,
                                   fromImageStatGenOuter = TRUE,
                                   numOfImgs = numOfImgs,
                                   otherPlotDat = primaryImageRaw,
                                   returnMat = FALSE)
  names(secondImageDat[[1]]) <- paste0(names(secondImageDat[[1]]), "_reference")
  names(primaryImageDat) <- paste0(names(primaryImageDat), "_focus")
  c(secondImageDat[[1]], primaryImageDat)
}