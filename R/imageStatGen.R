#' Generation of cell-based assay scores
#' 
#' This function is the central umbrella function that takes a list of image 
#' files or image file directories and generates a list of statistics for each
#' of the files. Needed input is the directory, the frame number and the nucleus
#' size cutoff. The latter can be found using the sizeCutoffGen function, which
#' is the only other function that the user needs to be concerned about. 
#' In this function, we start by defining a background threshold. This should be
#' individual, as the frames can have considerably different background 
#' intensities. We here take advantage of the auto_thresh function in the 
#' autothresholdr package, that uses the same background thresholding algorithms
#' as imageJ. See: https://imagej.net/plugins/auto-threshold
#' @seealso \code{\link[autothresholdr]{auto_thresh}}
#' @importFrom png readPNG writePNG
#' @importFrom tiff readTIFF
#' @importFrom RBioFormats read.image
#' @importFrom autothresholdr auto_thresh
#' @importFrom Matrix Matrix sparseMatrix summary
#' @importFrom stats sd median
#' @importFrom BiocParallel bplapply
#' @importFrom dbscan dbscan
#' @param imgDirs The directory containing the image files. Only nd2,
#' tiff and png formats are currently supported, and nd2 or non-normalised,
#' integer TIFF are clearly preferable.
#' There is an option here, which is mainly for example purposes, using the
#' provided dataset, and that is that a list of image files can be provided 
#' directly. These should then be in the form of a list of lists, each sublist 
#' being composed of matrices representing the colors. See example for usage. 
#' @param imgNames The names of the images. This argument is used both for
#' naming of the rows in the output matrix, but also for naming of the images
#' if such are generated using the diagnoImgs flag. 
#' @param frameNumIgG This identifies which of the frames in the files that 
#' contains the information about the autoantibody binding.
#' @param frameNumTrans Optionally, this identifies a frame in the files that 
#' contains the information about the transfection, i.e. often the GFP signal.
#' @param sizeCutoff The size of a typical cell. This is used to identify 
#' ring-shaped structures of a size equal to or larger than this cell threshold.
#' Smaller objects that are ring-shaped will be excluded with this cutoff, which
#' reduces noise considerably, but the threshold should not be set too high as
#' it then risks leading to the exclusion of true positive cells. This value
#' can be generated using the \code{\link{sizeCutoffGen}} function.
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. For reproducibility reasons, it might be clever to
#' run this multiple times, in this case, to bootstrap, or alternatively to
#' set a seed before starting. NB! If this command is used, then the images
#' will not be identical from round to round, and will not be comparable to the
#' output of other functions, as a randdom subset of the pictures are used. 
#' @param fixed Is the assay run with fixed or live cells? If the former, there
#' will be a second run of noise reduction. Not all fixed assays need this extra
#' precaution, but it is worth testing if the sesitivity is very poor. 
#' @param intensityCutoffIgG If an external intensity cutoff is used, e.g. if a 
#' large number of images from the same series are used, one can save
#' computational time by defining the background on the positive control. See
#' \code{\link{intensityCutoffGen}}. The two possible values here are a number
#' and "TRUE", in which each individual frame gets a separate intensity cutoff.
#' This latter option is suitable when it is unclear if all samples in a series
#' have been generated with the same values, or, more commonly, if RGB-
#' compliant files are used as input, as these are always normalised to a 
#' range from 0 to 1, and thus, the background in a negative frame will have 
#' considerably higher values than in a positive sample, simply because the 
#' rate from the lowest to the highest signal is more compressed in a case 
#' with noise only. 
#' @param intensityCutoffTrans Same as above but for the transfection control,
#' if present. 
#' @param diagnoImgs Should a images delimiting the islands that have been
#' selected be returned?
#' @param outDir The directory that the diagnoImages should be saved in. Only
#' used together with diagnoImgs = TRUE.
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @return A data frame with a number of statistics for the individual images.
#' Optionally, a directory containing images identifying the regions that have
#' been identified as cells of interest. 
#' @export imageStatGen
#' @examples
#' #Retrieve the example data
#' data(posImage)
#' data(noisyNegImage)
#' data(clearNegImage)
#' 
#' sizeCutoff <- sizeCutoffGen(imgDirs = clearNegImage, frameNum = 2)
#' intensityCutoff <- intensityCutoffGen(imgDirs = posImage, frameNum = 1)
#' 
#' #And now, for the final generation of results
#' allImgs <- list(posImage,noisyNegImage,clearNegImage)
#'                        
#' 
#' \dontrun{
#' result <- imageStatGen(imgDirs = allImgs,
#'                        imgNames =c("Pos", "Noisy_neg", "Clear_neg"),
#'                        frameNum = 1,
#'                        sizeCutoff = sizeCutoff,
#'                        intensityCutoff = intensityCutoff)
#'                        
#' #And here an example of how to use it with generation of images. 
#' #If the imgDirs is a list of nd2 directories, they can be converted to
#' #imgNames using this code: 
#' #imgNames <- gsub(".+/|.nd2", "", imgDirs)
#'
#' result <- imageStatGen(imgDirs = allImgs,
#'                        imgNames =c("Pos", "Noisy_neg", "Clear_neg"),
#'                        frameNum = 1,
#'                        sizeCutoff = sizeCutoff,
#'                        diagnoImgs = TRUE,
#'                        outDir = ".", 
#'                        intensityCutoff = intensityCutoff)
#' }
imageStatGen <- function(imgDirs, imgNames, frameNumIgG,
                         frameNumTrans = FALSE,
                         sizeCutoff, diagnoImgs = FALSE,
                         outDir = ".",
                         numPix = "All", fixed = FALSE, 
                         intensityCutoffIgG = TRUE, 
                         intensityCutoffTrans = TRUE,
                         numOfImgs = "All"){
  #First, we deal with the special case where the imgDir is not a directory
  #or a list of files, but an individual file, either containing only one
  #color matrix or a list of color matrices. 
  if(is.matrix(imgDirs)){
    imgDirs <- list(imgDirs)
  }
  #In the case where there was only one color matrix in the first place, or
  #where the file is made up of alist of color matrices, we go on to deal with
  #these here. 
  if(is.matrix(imgDirs[[1]])){
    imgDirs <- list(imgDirs)
  }
  #if(parallel){
  #  BPOPTIONS <- bpoptions()
  #} else {
  #  BPOPTIONS <- bpoptions(workers = 1)
  #}
  if(frameNumTrans == FALSE){
    resMat <- do.call("rbind", lapply(1:length(imgDirs), 
                                        imageStatGenInner, 
                                        imgDirs = imgDirs, 
                                        frameNum = frameNumIgG,
                                        sizeCutoff = sizeCutoff, 
                                        diagnoImgs = diagnoImgs, 
                                        outDir = outDir, 
                                        imgNames = imgNames,
                                        numPix = numPix, 
                                        fixed = fixed, 
                                        intensityCutoff = intensityCutoffIgG,
                                        numOfImgs = numOfImgs,
                                        returnMat = FALSE))
  } else {
    resMat <- do.call("rbind", lapply(1:length(imgDirs), imageStatGenOuter, 
                                        imgDirs = imgDirs, 
                                        frameNumIgG = frameNumIgG, 
                                        frameNumTrans = frameNumTrans,
                                        sizeCutoff = sizeCutoff, 
                                        diagnoImgs = diagnoImgs, 
                                        outDir = outDir, 
                                        imgNames = imgNames,
                                        numPix = numPix, 
                                        fixed = fixed, 
                                        intensityCutoffIgG = 
                                          intensityCutoffIgG,
                                        intensityCutoffTrans = 
                                          intensityCutoffTrans,
                                        numOfImgs = numOfImgs))
  }

  row.names(resMat) <- imgNames
  resMat
}