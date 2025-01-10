# Here, we calculate statistics on the island diameter: either the
# lowest or the mean diameter in the x-y directions. For each cell,
# that means the minimum of the number of columns and the number of
# rows that it distributes over.
islandDiameter <- function(islandPicture, statistic = "mean") {
    sm <- as.data.frame(Matrix::summary(islandPicture))
    colnames(sm) <- c("row", "column", "value")
    locRes <- unlist(lapply(split(sm, f = sm$value), function(x) {
        rowWidth <- max(x$row) - min(x$row) + 1
        colWidth <- max(x$col) - min(x$col) + 1
        if (statistic == "min") {
            min(c(rowWidth, colWidth))
        } else {
            mean(c(rowWidth, colWidth))
        }
    }))
}
