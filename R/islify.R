#' Generation of fluorescence-related data.
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
#' can be generated using the \code{\link{getSizeCutoff}} function.
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
#' \code{\link{getIntensityCutoff}}. The two possible values here are a number
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
#' 
#' 
#' 
#' See: \url{https://www.example.com} #This needs to be updated to 
#' 
#' 
#' 
#' 
#' 
#' @param intensityCutoffReference Same as above but for the control frame,
#' if present. 
#' @param fraction Within each 
#' @param reportIntensity Should the sum of the intensities for all the 
#' surviving islands be returned? Default is FALSE. 
#' @param diagnoImgs Should a images delimiting the islands that have been
#' selected be returned?
#' @param truncTo If diagnoImgs is TRUE, above which value should the data
#' be truncated? This argument has two possible inputs: a numeric value or
#' "max". Vectors of individual values, the same length as frameCols are also
#' accepted. The value option is often the most suitable, but requires some
#' knowledge of the range of values expected in the file. For this, the 
#' \code{\link{getQuantileIntensities}} function can be useful to provide
#' reasonable values. If this argument is set to a lower value than the true
#' max, the data will be truncated and a warning thrown.
#' @param outDir The directory that the diagnoImages should be saved in. Only
#' used together with diagnoImgs = TRUE.
#' @param numOfImgs If the provided files are nd2 format, they can contain
#' multiple files. In this case, this flag can be used to restrict the number of
#' used images.
#' @seealso \code{\link{getIntensityCutoff}}, \code{\link{getSizeCutoff}}
#' @return A data frame with statistics for the individual images.
#' Optionally, a directory containing images identifying the regions that have
#' been identified as cells of interest. 
#' @export islify
#' @examples
#' #Retrieve the example data
#' data(posImage)
#' data(negImage)
#' 
#' #First, establish the average nuclear size. Here, tje nuclei from one image
#' #is generally enough. 
#' sizeCutoff <- getSizeCutoff(imgDirs = list(negImage), frameNum = 3)
#' 
#' #We can run this algorithm in two ways: either only using the blue for size
#' #and red for the rest, or we can also incorporate the green color, saying 
#' #e.g.that we only are interested in IgG expression that colocalizes with GFP 
#' #expression. We start with the simpler case.
#' 
#' intensityCutoffRed <- getIntensityCutoff(imgDirs = list(negImage, posImage), 
#' frameNum = 1)
#' result <- islify(imgDirs = list(negImage, posImage),
#'                        imgNames =c("Neg", "Pos"),
#'                        frameNumFocus = 1,
#'                        diagnoImgs = FALSE,
#'                        sizeCutoff = sizeCutoff,
#'                        intensityCutoffFocus = 
#'                        intensityCutoffRed)
#'                        
#' #As can be noted above, diagnoImgs are set to FALSE, which means that the 
#' #output will be restricted to statistics only. If set to TRUE, then an image
#' #showing the selected areas in red on a black-and-white background are 
#' #generated. In that case, a directory, outDir, also needs to be specified. 
#' 
#' #Now comes the more complex useage, where the reference color is also 
#' #integrated. Here, to save time, only the negative image is used for 
#' #reference as both are expected to have similar GFP intensities. 
#' 
#' intensityCutoffGreen <- getIntensityCutoff(imgDirs = list(negImage), 
#' frameNum = 2)
#' 
#' resultWithGreen <- islify(imgDirs = list(negImage, posImage),
#'                        imgNames =c("Neg", "Pos"),
#'                        frameNumFocus = 1,
#'                        frameNumReference = 2,
#'                        diagnoImgs = FALSE,
#'                        sizeCutoff = sizeCutoff,
#'                        intensityCutoffFocus = 
#'                        intensityCutoffRed,
#'                        intensityCutoffReference = 
#'                        intensityCutoffGreen)
islify <- function(imgDirs, imgNames, frameNumFocus,
                         frameNumReference = FALSE,
                         sizeCutoff, diagnoImgs = FALSE,
                         truncTo = "max",
                         outDir = ".",
                         highNoise = FALSE, 
                         intensityCutoffFocus = TRUE, 
                         intensityCutoffReference = TRUE,
                         fraction = 0.5,
                         reportIntensity = FALSE,
                         numPix = "All", 
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
    resDf <- as.data.frame(do.call("rbind", lapply(1:length(imgDirs), 
                                      islifyInner, 
                                      imgDirs = imgDirs, 
                                      frameNum = frameNumFocus,
                                      sizeCutoff = sizeCutoff, 
                                      diagnoImgs = diagnoImgs, 
                                      truncTo = truncTo,
                                      outDir = outDir, 
                                      imgNames = imgNames,
                                      numPix = numPix, 
                                      highNoise = highNoise, 
                                      intensityCutoff = intensityCutoffFocus,
                                      fraction = fraction,
                                      reportIntensity = reportIntensity,
                                      numOfImgs = numOfImgs,
                                      returnMat = FALSE)))
    resDf$focFracTot <- resDf$totalPositivePixels/
    resDf$focFracTotDensAdjust <- 
      resDf$totalPositivePixels/resDf$islandDensity
    resDf <- resDf[,c(1,2,4)]
    colnames(resDf) <- c("intensityCutoff", "focFracTot",	
                         "focFracTotDensAdjust")
  } else {
    resDf <- as.data.frame(do.call("rbind", lapply(1:length(imgDirs), 
                                                   islifyOuter, 
                                                   imgDirs = imgDirs, 
                                                   frameNumFocus = 
                                                     frameNumFocus, 
                                                   frameNumReference = 
                                                     frameNumReference,
                                                   sizeCutoff = sizeCutoff, 
                                                   diagnoImgs = diagnoImgs, 
                                                   truncTo = truncTo,
                                                   outDir = outDir, 
                                                   imgNames = imgNames,
                                                   numPix = numPix, 
                                                   highNoise = highNoise, 
                                                   intensityCutoffFocus = 
                                                     intensityCutoffFocus,
                                                   fraction = fraction,
                                                   reportIntensity =
                                                     reportIntensity,
                                                   intensityCutoffReference = 
                                                     intensityCutoffReference,
                                                   numOfImgs = numOfImgs)))
  }

  row.names(resDf) <- imgNames
  resDf
}