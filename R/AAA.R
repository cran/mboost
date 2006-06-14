
### attach package `modeltools'

.onLoad <- function(lib, pkg) {
    if (!require("party"))
        stop("cannot load ", sQuote("party"))
    return(TRUE)
}
