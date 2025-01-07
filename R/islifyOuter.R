islifyOuter <- function(imgNum, imgDirs, frameNumFocus, frameNumReference,
                              sizeCutoff, 
                              diagnoImgs, outDir, imgNames, numPix, highNoise, 
                              intensityCutoffFocus, 
                              intensityCutoffReference,
                              fraction,
                              reportIntensity,
                              numOfImgs){
  #We start by analysing the reference picture
  secondImageDat <- islifyInner(imgNum = imgNum,
                                imgDirs = imgDirs, 
                                frameNum = frameNumReference,
                                sizeCutoff = sizeCutoff, 
                                diagnoImgs = FALSE, 
                                truncTo = FALSE,
                                numPix = numPix, 
                                highNoise = highNoise, 
                                intensityCutoff = intensityCutoffReference,
                                fraction = 1,
                                reportIntensity = FALSE,
                                fromIslifyOuter = "first",
                                numOfImgs = numOfImgs,
                                returnMat = TRUE)
  #Now, we locally import the focus image
  primaryImageRaw <- importFile(imgDirs[[imgNum]], frameNumFocus, numOfImgs)
  primaryImageFiltered <- primaryImageRaw
  primaryImageFiltered[which(secondImageDat[[2]] == 0)] <- 0
  primaryImageDat <- islifyInner(imgNum = 1,
                                   imgDirs = 
                                     list(abind(primaryImageFiltered, 
                                           matrix(0, 
                                                  nrow(primaryImageFiltered),
                                                  ncol(primaryImageFiltered)),
                                           along = 3)), 
                                   frameNum = 1,
                                   sizeCutoff = sizeCutoff, 
                                   diagnoImgs = diagnoImgs, 
                                   truncTo = truncTo,
                                   outDir = outDir, 
                                   imgNames = imgNames[[imgNum]],
                                   numPix = numPix, 
                                   highNoise = highNoise, 
                                   intensityCutoff = intensityCutoffFocus,
                                   fraction = fraction,
                                   reportIntensity = reportIntensity,
                                   fromIslifyOuter = "second",
                                   numOfImgs = numOfImgs,
                                   otherPlotDat = primaryImageRaw,
                                   returnMat = FALSE)
  names(secondImageDat[[1]]) <- paste0(names(secondImageDat[[1]]), "_reference")
  names(primaryImageDat) <- paste0(names(primaryImageDat), "_focus")
  c(secondImageDat[[1]], primaryImageDat)
}