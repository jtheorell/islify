#' Identification of an intensity cutoff for a whole experiment
#'
#' With this function, the otherwise integrated intensity cutoff is
#' externalised. This is only useful in situations where the original values
#' from the imaging have not been perturbed, so that the individual frames have
#' different max values, and are not all compressed to a range from 0 to 1. NB!
#' This requirement is not fulfilled if the files are exported as RGB compliant,
#' as the RGB format always ranges from 0 to 1.
#' @param imgDirs A vector or list of pathways, including filenames, to the
#' images to be analysed, e.g. "Raw_images/Positive_ctrl.nd2". Formats that are
#' currently supported are nd2, czi, tiff, png or lists of images, in the form
#' of three-dimensional arrays, where each layer in the third dimension
#' represents a color. nd2, czi or non-normalised, integer TIFF are clearly
#' preferable for memory and resolution purposes.
#' @param frameNum This identifies which of the frames in the file that
#' contains the information about the autoantibody binding.
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number
#' of used images.
#' @param numPix If the frames are very large, this can be used to reduce the
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting.
#' @param threshold_method The method used for thresholding. Available
#' alternatives are the same as for the
#' \code{\link[autothresholdr]{auto_thresh}} function. The default "Triangle"
#' method is in no way the only option, but it seems to perform reasonably
#' well under many circumstances.
#' @param ignore_white This is passed on to the
#' \code{\link[autothresholdr]{auto_thresh}} function. If a value, for example
#' inherited from getQuantileIntensities, then the values above this level will
#' not be considered when identifying the background threshold. Remedies some
#' of the variance between negative and positive control samples.
#' @seealso \code{\link[autothresholdr]{auto_thresh}}
#' @return An intensity cutoff value for the chosen color.
#' @examples
#' # Load example data and run the function:
#' data(negImage)
#' data(posImage)
#' getIntensityCutoff(imgDirs = list(negImage, posImage), frameNum = 1)
#' @export getIntensityCutoff
getIntensityCutoff <- function(imgDirs, frameNum,
                               numOfImgs = "All",
                               numPix = "All",
                               threshold_method = "Triangle",
                               ignore_white = FALSE) {
    # First, we deal with the special case where the imgDir is not a directory
    # or a list of files, but an individual file, either containing only one
    # color matrix or a list of color matrices.
    if (is.matrix(imgDirs)) {
        imgDirs <- list(imgDirs)
    }
    # In the case where there was only one color matrix in the first place, or
    # where the file is made up of a list of color matrices, we go on to deal 
    # with these here.
    if (is.matrix(imgDirs[[1]])) {
        imgDirs <- list(imgDirs)
    }
    median(unlist(lapply(imgDirs, function(x) {
        locFile <- importFile(x, frameNum, numOfImgs)
        if (is.numeric(numPix)) {
            locFile <- pixReduction(locFile, numPix)
        }
        # Now, we have to make a bit of a workaround again for the files that
        # are pre-normalised to the 0-1 range
        if (max(locFile) <= 1) {
            if(is.numeric(ignore_white) && ignore_white <= 1){
                ignore_white <- ignore_white * 1000
            }
            auto_thresh(round(locFile * 1000), threshold_method,
                ignore_white = ignore_white
            ) / 1000
        } else {
            auto_thresh(locFile,
                method = threshold_method,
                ignore_white = ignore_white
            )
        }
    })))
}
