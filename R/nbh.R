# generic function for nbh that supports object of class "integer" and "GRanges"
nbh <- function(x, ...) UseMethod("nbh")
