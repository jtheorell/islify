rescaleFrom01 <- function(locFile, intensityCutoff, printMessage = TRUE){
  if(printMessage){
    message("The max value in this file is 1. Therefore, it will be re-scaled",
            " and rounded, so that the max value is 1000, and all values are ",
            "integers.")
  }
  if(is.numeric(intensityCutoff)){
    if(intensityCutoff > 1){
      stop("The provided intensityCutoff is higher than the max value, so all",
           " datapoints will be removed. Change this and try again.")
    } else {
      intensityCutoff <- intensityCutoff*1000
    }
  }
  #Here, we need to change the values into 
  locFile <- round(locFile*1000)
  list(locFile, intensityCutoff)
}