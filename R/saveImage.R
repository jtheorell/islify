#' Saving individually coloured png images from raw image input
#'
#' This function simply generates pngs from raw image files, such as .nd2.
#'
#' @param imgDirs The directory containing the image files. nd2, czi,
#' tiff and png formats or a list of three-dimensional arrays, where each
#' layer in the third dimension represents a color are currently supported.
#' nd2, czi or non-normalised, integer TIFF are clearly preferable.
#' @param imgNames The names of the images. This argument is used both for
#' naming of the rows in the output matrix, but also for naming of the images
#' if such are generated using the diagnoImgs flag.
#' @param frameNums This identifies which of the frames in the files that
#' should be plotted. Defaults to "All", which then also generates a merge.
#' @param frameCols Which color does the frame/frames in question have?
#' Strangely, in some file formats such as nd2, it varies - not all are
#' "R-G-B", but "G-R-B" also occurs. Can take the values "R"/"Red", 
#' "G"/"Green", and "B"/"Blue". Should be a vector the same length as 
#' frameNums. If nothing is provided, it defaults to only saving the merge.
#' @param truncTo Above which value should the data be truncated? This argument
#' has two possible inputs: a numeric value or "max". Vectors of
#' individual values, the same length as frameCols are also accepted. The value
#' option is often the most suitable, but requires some knowledge of the range
#' of values expected in the file. For this, the
#' \code{\link{getQuantileIntensities}} function can be useful to provide
#' insights to provide reasonable values. If this value is set to a lower value
#' than the true max, the data will be truncated and a warning thrown.
#' @param outDir The directory that the images should be saved in. Subfolders
#' are created within this directory for each color.
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number
#' of used images.
#' @param numPix If the frames are very large, this can be used to reduce the
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. NB! If this command is used, then the images
#' will not be identical from round to round, and will not be comparable to the
#' output of other functions, as a random subset of the pictures are used.
#' @param saveImages Should images be saved? This of course seems like a rather
#' pointless command for a function meant to save pictures, but it is there
#' for development reasons. Don't use it!
#' @seealso \code{\link{getQuantileIntensities}}
#' @return Named png images that are scaled to the max allowed value in the
#' raw image file and coloured according to frameCol.
#' @examples
#' # Retrieve the example data
#' data(posImage)
#' data(negImage)
#'
#' # To use this function, the getQuantileIntensities function is very helpful,
#' # as it can be used to define reasonable truncTo values. It is of course
#' # important to use a strongly stained sample in this case, to make sure that
#' # the positive samples in the series are not overly truncated, reducing the
#' # visual signal-to-noise ratio.
#' posQuants <- getQuantileIntensities(list(posImage), quantiles = 0.99)
#' # We will use the 99th percentile
#' posQuants
#' #    Percent_99
#' # 1  0.3681333
#' # 2  0.2476340
#' # 3  0.3794010
#'
#' # And now to the function. NB! In this case, we are using the last
#' # "saveImages = FALSE" flag, which means that nothing is produced! Remove it
#' # to use the function.
#' saveImage(list(posImage, negImage), c("pos", "neg"),
#'     frameNums = "All",
#'     frameCols = c("R", "G", "B"), truncTo = posQuants[[1]][, "Percent_99"],
#'     outDir = ".", saveImages = FALSE
#' )
#'
#' @export saveImage
saveImage <- function(imgDirs, imgNames, frameNums = "All", frameCols,
                      truncTo = "max", outDir = ".", numOfImgs = "All",
                      numPix = "All", saveImages = TRUE) {
    if (saveImages) {
        if (frameNums == "All") {
            for (i in c("Merged_colors", frameCols)) {
                dir.create(paste0(outDir, "/", i))
            }
        }
        lapply(seq_along(imgDirs), function(x) {
            imgDir <- imgDirs[[x]]
            imgName <- imgNames[[x]]
            # We start by the default, i.e. that we want to plot everything.
            if (frameNums == "All") {
                locFile <- importFile(imgDir,
                    frameNum = "All",
                    numOfImgs = numOfImgs
                )
                # Here, if max is 1, all other normalisation
                # information will be discarded with a warning.
                if (max(locFile) == 1 && (truncTo[1] > 1)) {
                    message(
                        "This file has a max value of 1, and therefore this",
                        " will be used for truncation purposes instead of the", 
                        " truncTo value. "
                    )
                    truncTo <- 1
                }

                if (is.numeric(numPix)) {
                    locFile <- pixReduction(locFile, numPix)
                }
                if (length(truncTo) == 1) {
                    if (truncTo == "max") {
                        truncTo <- vapply(
                            seq_along(frameCols),
                            function(y) max(locFile[, , y]), 1
                        )
                    } else {
                        truncTo <- rep(truncTo, length(frameCols))
                    }
                }
                if (length(truncTo) != length(frameCols)) {
                    stop(
                        "Mismatch between the length of frameCols and truncTo",
                        ". Change this. "
                    )
                }
                # First, we make the plots for individual colors
                for (i in seq_along(frameCols)) {
                    frameCol <- frameCols[i]
                    saveImageInner(
                        locImg = locFile[, , i],
                        frameCol = frameCol,
                        truncToInner = truncTo[i],
                        imgName = imgName,
                        outDir = outDir
                    )
                }
                # Then the common plot.
                # Here, we normalise to the range of the possible max value.
                locImgNorm <- locFile
                for (i in seq_along(frameCols)) {
                    locImgNorm[, , i] <- locFile[, , i] / truncTo[i]
                }

                frameNums <- vapply(frameCols, colToNum, 1)
                if (max(locImgNorm) > 1) {
                    locImgNorm[which(locImgNorm > 1)] <- 1
                }
                # And now, we create a full rgb array.
                locArray <- array(0, dim = c(
                    nrow(locImgNorm),
                    ncol(locImgNorm),
                    3
                ))
                for (i in seq_along(frameCols)) {
                    locArray[, , frameNums[i]] <- locImgNorm[, , i]
                }

                writePNG(locArray, paste0(
                    outDir, "/Merged_colors/",
                    imgName, "_merge.png"
                ))
            } else {
                locFile <- importFile(imgDir,
                    frameNum = frameNums,
                    numOfImgs = numOfImgs
                )
                saveImageInner(
                    locImg = locFile,
                    frameCol = frameCol,
                    truncToInner = truncTo,
                    imgName = imgName,
                    outDir = outDir
                )
            }
        })
    }
}
