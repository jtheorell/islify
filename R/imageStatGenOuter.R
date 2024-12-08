imageStatGenOuter <- function(imgNum, imgDirs, frameNumIgG, frameNumTrans,
                              sizeCutoff, 
                              diagnoImgs, outDir, imgNames, numPix, fixed, 
                              intensityCutoffIgG, 
                              intensityCutoffTrans, 
                              numOfImgs){
  #We start by analysing the transfection picture
  transImageDat <- imageStatGenInner(imgNum = imgNum,
                                     imgDirs = imgDirs, 
                                     frameNum = frameNumTrans,
                                     sizeCutoff = sizeCutoff, 
                                     diagnoImgs = FALSE, 
                                     numPix = numPix, 
                                     fixed = fixed, 
                                     intensityCutoff = intensityCutoffTrans,
                                     numOfImgs = numOfImgs,
                                     returnMat = TRUE)
  #Now, we locally import the IgG image
  igGImageRaw <- importFile(imgDirs[[imgNum]], frameNumIgG, numOfImgs)
  igGImageFiltered <- igGImageRaw
  igGImageFiltered[which(transImageDat[[2]] == 0)] <- 0
  igGImageDat <- imageStatGenInner(imgNum = 1,
                                   imgDirs = list(list(igGImageFiltered)), 
                                   frameNum = 1,
                                   sizeCutoff = sizeCutoff, 
                                   diagnoImgs = diagnoImgs, 
                                   outDir = outDir, 
                                   imgNames = imgNames[[imgNum]],
                                   numPix = numPix, 
                                   fixed = fixed, 
                                   intensityCutoff = intensityCutoffIgG,
                                   numOfImgs = numOfImgs,
                                   otherPlotDat = igGImageRaw,
                                   returnMat = FALSE)
  allImageDat <- c(transImageDat[[1]], igGImageDat)
  names(allImageDat) <- paste0(names(allImageDat),
                               rep(c("_transfection", "_IgG"), 
                                   each = length(igGImageDat)))
  allImageDat
}