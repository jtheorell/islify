#' Generation of fluorescence-related scores.
#' 
#' This function is the central umbrella function that takes a list of image 
#' files or image file directories and generates a list of statistics for each
#' of the files. 
#' @seealso \code{\link[autothresholdr]{auto_thresh}}
#' @importFrom png readPNG writePNG
#' @importFrom tiff readTIFF
#' @importFrom RBioFormats read.image
#' @importFrom autothresholdr auto_thresh
#' @importFrom Matrix Matrix sparseMatrix summary
#' @importFrom stats sd median quantile
#' @importFrom dbscan dbscan
#' @importFrom abind abind
#' @importFrom methods new 
#' @param imgDirs A vector or list of pathways, including filenames, to the
#' images to be analysed, e.g. "Raw_images/Positive_ctrl.nd2". Formats that are 
#' currently supported are nd2, czi, tiff, png or lists of images, in the form 
#' of three-dimensional arrays, where each layer in the third dimension 
#' represents a color. nd2, czi or non-normalised, integer TIFF are clearly 
#' preferable for memory and resolution purposes. 
#' @param imgNames The names of the images. This argument is used both for
#' naming of the rows in the output matrix, but also for naming of the images
#' if such are generated using the diagnoImgs flag. 
#' @param frameNumFocus This identifies which of the frames in the files that 
#' contains the information that should be measured primarily, such as IgG
#' binding in the case of autoantibody screening. 
#' @param frameNumReference Optionally, this identifies a frame in the files 
#' that can be used to focus the analysis to certain regions. It could be
#' that only areas with GFP expression can contain meaningful expression of 
#' the marker of primary interest, and then only these areas will be considered.
#' @param sizeCutoff The size of a typical cell nucleus. This is used to 
#' identify structures of a size equal to or larger than this cell threshold.
#' Smaller objects that are ring-shaped will be excluded with this cutoff, which
#' reduces noise considerably, but the threshold should not be set too high as
#' it then risks leading to the exclusion of true positive cells. This value
#' can be generated using the \code{\link{sizeCutoffGen}} function.
#' @param numPix If the frames are very large, this can be used to reduce the 
#' computational burden. NB! If this command is used, then the images
#' will not be identical from round to round, and will not be comparable to the
#' output of other functions, as a random subset of the pictures are used. 
#' @param highNoise Is the assay especially noisy, with very high and individual
#' background? If so, a second run of noise reduction can be used to increase
#' signal-to-noise ratio. This is often but not always needed for fixed 
#' cell-based assays, e.g.
#' @param intensityCutoffFocus If an external intensity cutoff is used, e.g. if a 
#' large number of images from the same series are used, one can save
#' computational time by defining the background on the positive control. See
#' \code{\link{intensityCutoffGen}}. The two possible values here are a number
#' and "TRUE", in which each individual frame gets a separate intensity cutoff.
#' This latter option is suitable when it is unclear if all samples in a series
#' have been generated with the same imaging settings, or, more commonly, if RGB-
#' compliant files are used as input, as these are always normalised to a 
#' range from 0 to 1, and thus, the background in a negative frame will have 
#' considerably higher values than in a positive sample, simply because the 
#' rate from the lowest to the highest signal is more compressed in a case 
#' with noise only. If individual values are used, thresholding is performed 
#' with the auto_threshold function from the autothresholdr package, that uses
#' the same thresholding algorithms as imageJ. Specifically, it is the "tri"
#' option that is used. 
#' See: \url{https://www.example.com}
#' @param intensityCutoffReference Same as above but for the control frame,
#' if present. 
#' @param reportIntensity Should the sum of the intensities for all the 
#' surviving islands be returned? Default is FALSE. 
#' @param nuclearToCellQuotient The quotient of the median nuclear width divided
#' by the median cell width. Set to 0.66, as this is the quotient for HEK cells
#' according to the literature. This is used to calculate how many cells that
#' are clumped together in large, confluent multicell objects in the pictures. 
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
#' @export getImageStats
#' @examples
#' #Retrieve the example data
#' data(posImage)
#' data(negImage)
#' 
#' #First, establish the average nuclear size. Here, tje nuclei from one image
#' #is generally enough. 
#' sizeCutoff <- getSizeCutoff(imgDirs = list(negImage), frameNum = 3)
#' 
#' #Now, we start with the optimal scenario, where there, in addition to the 
#' #blue nuclear dye and the e.g. red patient IgG dye, is also a transfection 
#' #control dye, such as GFP. In this case, the algorithm will filter out all 
#' #object with a smaller diameter (in straig x-x, y-y lines) in the 
#' #transfection control dye than the nuclear width. Then it identifies objects
#' #that are stained with e.g. the IgG dye and that are wider and highe than
#' #the standard nuclear width.
#' 
#' 
#' we can run this in two ways: either only using the blue for size and 
#' #red for the rest, or we can also incorporate the green color, saying e.g.
#' #that we only are 
#' 
#' 
#' 
#' 
#' intensityCutoffRed <- getIntensityCutoff(imgDirs = list(negImage, posImage), 
#' frameNum = 1)
#' 
#' #And now, for the final generation of results
#' result <- getImageStats(imgDirs = list(negImage, posImage),
#'                        imgNames =c("Neg", "Pos"),
#'                        frameNum = 1,
#'                        diagnoImgs = FALSE,
#'                        sizeCutoff = sizeCutoff,
#'                        intensityCutoff = intensityCutoff)

getImageStats <- function(imgDirs, imgNames, frameNumFocus,
                         frameNumReference = FALSE,
                         sizeCutoff, diagnoImgs = FALSE,
                         outDir = ".",
                         numPix = "All", highNoise = FALSE, 
                         intensityCutoffFocus = TRUE, 
                         intensityCutoffReference = TRUE,
                         reportIntensity = FALSE,
                         nuclearToCellQuotient = 0.66,
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

  if(frameNumReference == FALSE){
    resMat <- do.call("rbind", lapply(1:length(imgDirs), 
                                      imageStatGenInner, 
                                      imgDirs = imgDirs, 
                                      frameNum = frameNumFocus,
                                      sizeCutoff = sizeCutoff, 
                                      diagnoImgs = diagnoImgs, 
                                      outDir = outDir, 
                                      imgNames = imgNames,
                                      numPix = numPix, 
                                      highNoise = highNoise, 
                                      intensityCutoff = intensityCutoffFocus,
                                      reportIntensity = reportIntensity,
                                      nuclearToCellQuotient = 
                                        nuclearToCellQuotient,
                                      numOfImgs = numOfImgs,
                                      returnMat = FALSE))
  } else {
    resMat <- do.call("rbind", lapply(1:length(imgDirs), imageStatGenOuter, 
                                      imgDirs = imgDirs, 
                                      frameNumFocus = frameNumFocus, 
                                      frameNumReference = frameNumReference,
                                      sizeCutoff = sizeCutoff, 
                                      diagnoImgs = diagnoImgs, 
                                      outDir = outDir, 
                                      imgNames = imgNames,
                                      numPix = numPix, 
                                      highNoise = highNoise, 
                                      intensityCutoffFocus = 
                                        intensityCutoffFocus,
                                      reportIntensity =
                                        reportIntensity,
                                      nuclearToCellQuotient = 
                                        nuclearToCellQuotient,
                                      intensityCutoffReference = 
                                        intensityCutoffReference,
                                      numOfImgs = numOfImgs))
  }

  row.names(resMat) <- imgNames
  as.data.frame(resMat)
}