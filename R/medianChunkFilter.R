medianChunkFilter <- function(locFileRaggedIsles, locFile, threshold_method,
                              ignore_white){
    origDims <- dim(locFile)
    
    #Now, which side is longer? 
    fileRows <- nrow(locFile)
    fileCols <- ncol(locFile)
    rowColRatio <- fileRows/fileCols
    
    if(rowColRatio < 1){
        nChunksRow <- 4
        nChunksCol <- round(4/rowColRatio)
    } else {
        nChunksRow <- round(4*rowColRatio)
        nChunksCol <- 4
    }
    rowNumList <- split(seq(1, fileRows), 
                        ceiling(seq(1, fileRows)/(fileRows/nChunksRow)))
    colNumList <- split(seq(1, fileCols), 
                        ceiling(seq(1, fileCols)/(fileCols/nChunksCol)))
    #Now, we calculate the background for each of the squares
    backgroundMat <- matrix(NA, nChunksRow, nChunksCol)
    for(i in seq_along(rowNumList)){
        for(j in seq_along(colNumList)){
            rows <- rowNumList[[i]]
            cols <- colNumList[[j]]
            backgroundMat[i,j] <- auto_thresh(locFile[rows,cols], 
                                         method = "Triangle",
                                         ignore_white = ignore_white)
        }
    }
    #And here the background for all, as the individual backgrounds tend 
    #to be variable and in some cases very low, dragging down the median
    #to impractical levels. 
    backgroundAll <- auto_thresh(locFile, 
                                 method = "Triangle",
                                 ignore_white = ignore_white)
    backgroundMad <- mad(backgroundMat)
    locFileDark <- locFileRaggedIsles
    for(i in seq_along(rowNumList)){
        for(j in seq_along(colNumList)){
            rows <- rowNumList[[i]]
            cols <- colNumList[[j]]
            locBackground <- backgroundMat[i,j]
            if(locBackground > backgroundAll+(1*backgroundMad)){
                locFileDark[rows,cols] <- 0
            }
        }
    }
    
    bigIslesPixels <- islandPixels(locFileRaggedIsles)
    darkIslesPixels <- islandPixels(locFileDark)
    if(length(bigIslesPixels) > length(darkIslesPixels)){
        overlappingIsles <- which(names(bigIslesPixels) %in%
                                      names(darkIslesPixels))
        missingIsles <- 
            names(bigIslesPixels)[-overlappingIsles]
        bigIslePixelsRed <- bigIslesPixels[overlappingIsles]
        smallerIsles <- names(bigIslesPixels)[which(unlist(bigIslePixelsRed) > 
                                                 unlist(darkIslesPixels))]
        darkIsles <- c(missingIsles, smallerIsles)
    } else {
        darkIsles <- names(bigIslesPixels)[which(unlist(bigIslesPixels) > 
                                              unlist(darkIslesPixels))]
    }
    
    if(length(darkIsles) > 0){
        #Now,  we create the new, reduced file. 
        sm <- as.data.frame(summary(locFileRaggedIsles))
        colnames(sm) <- c("row", "column", "value")
        sm$value[which(sm$value %in% as.numeric(names(darkIsles)))] <- 0
        # Here, we remove all new zeros.
        smSmall <- sm[which(sm$value > 0), ]
        locSM <- sparseMatrix(
            i = smSmall$row,
            j = smSmall$column,
            x = smSmall$value,
            dims = origDims
        )
    } else {
        locFileRaggedIsles
    }
}