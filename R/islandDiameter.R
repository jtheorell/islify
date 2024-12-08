#Now, we calculate the average diameter. For each cell, that means the average
#of the number of columns and the number of rows that it distributes over.
islandDiameter <- function(islandPicture){
    sm <- as.data.frame(Matrix::summary(islandPicture))
    colnames(sm) <- c("row", "column", "value")
    locRes <- unlist(lapply(split(sm, f = sm$value), function(x){
      rowWidth <- max(x$row)-min(x$row)+1
      colWidth <- max(x$col)-min(x$col)+1
      mean(c(rowWidth, colWidth))
    }))
}
