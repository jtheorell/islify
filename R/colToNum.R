colToNum <- function(frameCol) {
    switch(frameCol,
        "R" = 1,
        "Red" = 1,
        "red" = 1,
        "G" = 2,
        "Green" = 2,
        "green" = 2,
        "B" = 3,
        "Blue" = 3,
        "blue" = 3
    )
}
