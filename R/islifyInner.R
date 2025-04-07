islifyInner <- function(imgNum, imgDirs, frameNum, 
                        frameNumNuclei = FALSE, sizeCutoff,
                        diagnoImgs, truncTo, outDir, imgNames, numPix,
                        highNoise, intensityCutoff, 
                        intensityCutoffNuclei = FALSE,
                        threshold_method,
                        ringFrac, flatFrac, truncLim = "max", numOfImgs,
                        reportIntensity, ignore_white,
                        otherPlotDat = FALSE, fromIslifyOuter = FALSE) {
    imgDir <- imgDirs[[imgNum]]
    
    locFile <- importFile(imgDir, frameNum, numOfImgs,
        fromIslifyOuter = fromIslifyOuter
    )
    # In cases where the data already has been compressed, we need to
    # integerize it for the functions below to work.
    if (max(locFile) <= 1) {
        if (fromIslifyOuter == "first") {
            printMessage <- FALSE
        } else {
            printMessage <- TRUE
        }
        rescaledList <- rescaleFrom01(
            locFile,
            intensityCutoff,
            printMessage
        )
        locFile <- rescaledList[[1]]
        intensityCutoff <- rescaledList[[2]]
        if (is.numeric(truncTo) && truncTo <= 1) {
            truncTo <- truncTo * 1000
        }
        if(is.numeric(ignore_white) && ignore_white <= 1){
            ignore_white <- ignore_white*1000
        }
    }
    
    # Here, we restrict the number of pixels, in applicable cases, to a
    # randomly located subset of the frame
    if (is.numeric(numPix)) {
        locFile <- pixReduction(locFile, numPix)
    }

    locFileClean <- locFile

    if (is.logical(intensityCutoff)){
        if(intensityCutoff){
            intensityCutoff <- as.numeric(auto_thresh(locFile, 
                                                      threshold_method, 
                                                      ignore_white))
        } else {
            intensityCutoff <- min(locFile)
        }
    } 
    
    if(inherits(frameNumNuclei, what = "integer")){
        locFileNucleus <- importFile(imgDir, frameNumNuclei, numOfImgs,
                                     fromIslifyOuter = fromIslifyOuter
        )
        if(is.logical(intensityCutoffNuclei) && intensityCutoffNuclei == TRUE){
            intensityCutoffNuclei <- as.numeric(auto_thresh(locFileNucleus, 
                                                 threshold_method, 
                                                 ignore_white = ignore_white))
        }
        numOfNucleatedPixels <- 
            length(which(as.vector(locFileNucleus) > intensityCutoffNuclei))
    }

    #Here we define the truncation limit. 
    if(truncLim == "max"){
        truncLim <-  max(locFileClean)
    }
    
    locFileClean[which(locFile <= intensityCutoff)] <- 0
    
    if (length(which(locFileClean != 0)) / length(locFileClean) > 0.9) {
        warning(
            "More than 90% of the pixels are positive with this threshold. ",
            "This is unlikely to be correct even in a ",
            "positive sample or control. "
        )
    }
    # Now, we are going to introduce an optional second filtering step, where
    # we disregard the information from the areas devoid of cells, and also
    # remove the noise from the non-important cells.
    if (highNoise) {
        intensityCutoff2 <- as.numeric(auto_thresh(round(locFileClean * 100),
            threshold_method,
            ignore_black = TRUE,
            ignore_white = ignore_white
        ) / 100)
        locFileClean[which(locFile <= intensityCutoff2)] <- 0
    }

    
    resList <- list()
    # We include information about the first filter threshold
    resList$intensityCutoff <- intensityCutoff
    
    #Before going anywhere else, we will also truncate the values above the
    #truncation limit. 
    locFileClean[which(locFileClean > truncLim)] <- truncLim
    
    locFile01 <- locFileClean

    # In the case that there is nothing but noise below background, a negative
    # result is thrown back here
    if (length(which(locFile01 != 0)) == 0
        ) {
        resList$fractionOfAll_focus <- 0
        if(inherits(frameNumNuclei, what = "integer")){
            resList$fractionOfNuclei <- 0
        }
        if (reportIntensity) {
            resList$intensitySum <- 0
        }
    } else {
        locFile01[which(locFile01 != 0)] <- 1
        locFileIslified <- islandNaming(locFile01)

        # Importantly, we also need to remove the now zero points from the
        # locFileClean, as some points are excluded by the dbscan function.
        if (length(which(locFileIslified == 0)) > 0) {
            locFileClean[which(locFileIslified == 0)] <- 0
        }
        locFileIslifiedSparse <- as(locFileIslified, "dgCMatrix")

        diameters <- islandDiameter(locFileIslifiedSparse, statistic = "min")
        if (any(diameters > sizeCutoff)) {
            locFileBigIsles <- islandRemoval(
                locFileIslifiedSparse, diameters,
                sizeCutoff
            )
            # Here come a few steps that are not relevant if we are running the
            # outer function for a co-localised protein.
            if (fromIslifyOuter != "first") {
                locFileClean[which(as.matrix(locFileBigIsles) == 0)] <- 0
                locHighIsles <- islandPeaks(locFileBigIsles,
                    intensityVec =
                        as(locFileClean, "dgCMatrix")@x,
                    ringFrac = ringFrac,
                    flatFrac = flatFrac
                )
                diametersHigh <- islandDiameter(locHighIsles,
                    statistic = "min"
                )

                # Now, we are going to create an interesting criterion: only
                # islands that retain their diameter even after all values
                # below the fraction of the 99th percentile are kept.
                locFileRaggedIsles <- islandRemoval(
                    locHighIsles,
                    diametersHigh,
                    sizeCutoff
                )
                
                #And now a further criterion. Given that clumps of cells at
                #times give rise to considerable extra background, which in
                #all other ways "seems" like true signal to the software, 
                #we here remove any island that is part of a large chunk of
                #the picture where the median is higher than the median+3MAD
                #for the whole picture. THis only applies if the picture is
                #larger in one direction than 2000 pixels. 
                
                if(max(dim(locFileRaggedIsles)) > 2000){
                    locFileBigIsles <- medianChunkFilter(locFileRaggedIsles, 
                                                         locFile, 
                                                         threshold_method, 
                                                         ignore_white)
                } else {
                    locFileBigIsles <- locFileRaggedIsles
                }
                
                if (reportIntensity) {
                    # Here, we add information about the sum of
                    # intensity in the big islands
                    locFileBigIslesMat <- as.matrix(locFileBigIsles)
                    locFileClean[which(locFileBigIslesMat == 0)] <- 0
                    resList$intensitySum <- sum(locFileClean)
                }
            }
            # Now, we calculate the number of pixels per island
            bigIslandPixels <- islandPixels(locFileBigIsles)

            # Here, we export the number of positive pixels as a fraction of
            # the total number of pixels in the frame.
            resList$fractionOfAll_focus <-
                sum(unlist(bigIslandPixels)) /
                    (dim(locFileBigIsles)[1] * dim(locFileBigIsles)[2])
            
            #And if we have the nucleus information, this is exported here. 
            if(inherits(frameNumNuclei, what = "integer")){
                resList$fractionOfNuclei <- 
                    sum(unlist(bigIslandPixels))/numOfNucleatedPixels
            }
                
        } else {
            resList$fractionOfAll_focus <- 0
            if(inherits(frameNumNuclei, what = "integer")){
                resList$fractionOfNuclei <- 0
            }
            if (reportIntensity) {
                resList$intensitySum <- 0
            }
        }

        if (resList$fractionOfAll_focus > 0) {
            bigClean <- as.matrix(locFileClean)
            locFileBigIslesMat <- as.matrix(locFileBigIsles)
            bigClean[which(locFileBigIslesMat == 0)] <- 0
        } else {
            bigClean <- matrix(0, nrow(locFileClean), ncol(locFileClean))
        }
    }
    # Now, we include the second threshold
    if (highNoise) {
        resList$secondIntensityCutoff <- intensityCutoff2
    }
    # At this stage, we can, if it is prompted, throw back a picture of the
    # selected islands in the context of the full file.
    if (diagnoImgs) {
        if (is.logical(otherPlotDat) == FALSE && is(otherPlotDat, "matrix")) {
            if (max(otherPlotDat) <= 1) {
                otherPlotDat <- round(otherPlotDat * 1000)
            }
            locFile <- otherPlotDat
        }
        plotFile <- array(
            matrix(0, nrow(locFile), ncol(locFile)),
            c(nrow(locFile), ncol(locFile), 3)
        )
        # In the case where another file has been added that should be the
        # background plot, i.e. when we are using the imageStatGenOuter, we
        # add this here
        plotDat <- locFile

        if (truncTo == "max") {
            truncVal <- max(plotDat)
        } else {
            truncVal <- truncTo
            plotDat[which(plotDat > truncVal)] <- truncVal
        }
        locFile01 <- (plotDat - min(plotDat)) / (truncVal - min(plotDat))
        # And here, we invert the colors.
        locFile10 <- abs(locFile01 - 1)
        plotFile[, , 1] <- plotFile[, , 2] <- plotFile[, , 3] <- locFile10

        # And now to the magic: we will remove all values from the green and
        # blue frames, if there are any selected cells. Else, this step is
        # excluded.
        # We start by not fully flattening, but still decreasing the
        # range of the positive values.
        locFile0.4 <- (locFile01 * 0.6) + 0.4
        if (resList$fractionOfAll_focus > 0) {
            plotFile[, , 1][which(bigClean != 0)] <- 
                locFile0.4[which(bigClean != 0)]
            plotFile[, , 2][which(bigClean != 0)] <- 0
            plotFile[, , 3] <- plotFile[, , 2]
        } else {
            plotFile <- plotFile
        }
        writePNG(plotFile, paste0(outDir, "/", imgNames[imgNum], ".png"))
    }
    if (fromIslifyOuter == "first") {
        list(unlist(resList), locFileClean)
    } else {
        unlist(resList)
    }
}
