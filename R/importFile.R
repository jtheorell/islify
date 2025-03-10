importFile <- function(imgDir, frameNum, numOfImgs,
                       fromIslifyOuter = FALSE) {
    if (is.array(imgDir)) {
        # Here we send a message only in the cases where the reason
        # for this aberrance is not coming from within the software itself.
        if (fromIslifyOuter == "second") {
            message(
                "This object is interpreted as an image",
                " file and not a directory."
            )
        }
        if (is.numeric(frameNum)) {
            locFile <- imgDir[, , frameNum]
        } else if (frameNum == "All") {
            locFile <- imgDir
        }
    } else if (grepl("\\.nd2", imgDir)) {
        message(imgDir)
        # For this, we need to be clever, as many files can be saved in the
        # same object. For now, we will allow all to go through.
        locFileSet <- read.image(imgDir, normalize = FALSE)
        if (numOfImgs == "All") {
            numOfImgs <- length(locFileSet@.Data)
        }

        if (is.numeric(frameNum)) {
            locFile <- do.call("cbind", lapply(seq(1, numOfImgs), function(x) {
                locFileSet[[x]]@.Data[, , frameNum]
            }))
        } else if (frameNum == "All") {
            locFileList <- lapply(seq(1, numOfImgs), function(x) {
                locFileSet[[x]]@.Data
            })
            locFile <- abind(locFileList, along = 2)
        } else {
            stop("Cannot interpret frameNum. Change input to All or a value.")
        }
    } else if (grepl(".png", imgDir)) {
        message(imgDir)
        locFile <- readPNG(imgDir)
        if (is.numeric(frameNum)) {
            locFile <- locFile[, , frameNum]
        } else if (frameNum == "All") {
            locFile
        } else {
            stop("Cannot interpret frameNum. Change input to All or a value.")
        }
    } else if (grepl(".tif", imgDir)) {
        message(imgDir)
        locFile <- readTIFF(imgDir, as.is = TRUE, all = FALSE)
        if (is.numeric(frameNum)) {
            locFile <- locFile[, , frameNum]
        } else if (frameNum == "All") {
            locFile
        } else {
            stop("Cannot interpret frameNum. Change input to All or a value.")
        }
    } else if (grepl("\\.czi", imgDir)) {
        message(imgDir)
        locFile <- read.image(imgDir, normalize = FALSE)@.Data
        # This file format can contain essentially anything. For that reason,
        # we are conservative here, and will only allow an n-dimensional array
        # to go through.
        if (is.array(locFile) || is.matrix(locFile)) {
            if (is.numeric(frameNum)) {
                locFile <- locFile[, , frameNum]
            } else if (frameNum == "All") {
                locFile
            } else {
                stop("Cannot interpret frameNum.",
                " Change input to All or a value.")
            }
        }
    } else {
        stop(
            "This file format is not supported yet. Write to the ",
            "maintainer or change the format"
        )
    }
    locFile
}
