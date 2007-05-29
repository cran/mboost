
library("utils")
library("tools")
Rnw <- list.files(pattern = "Rnw")
sapply(Rnw, function(f) Sweave(f))
texi2dvi("spatial_boosting.tex", pdf = TRUE, clean = FALSE)
