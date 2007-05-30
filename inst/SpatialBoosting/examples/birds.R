
library("mboost")
attach(asNamespace("mboost"))

source("../R/spaboost.R")
library("mgcv")

load("../data/birds.Rda")

fm <- as.formula(paste("SG4 ~ ", paste(indep, collapse = " + ")))
sp <- spaboost(fm, data = birds, spatial = coord, family = Poisson())

indx <- which(sp$ensemble == 33)
spat <- rowSums(sapply(indx, function(i) 
    predict(sp$ensembless[[i]]) * sp$control$nu))

library("akima")
persp(interp(x = birds$x_gk, y = birds$y_gk, z = exp(spat)), 
    theta = -40, d = 200, ticktype = "det")

source("../R/gamboost2.R")
source("../R/base.R")

fm <- paste("SG4 ~ ", paste("mybs(", indep[1:4], ")", collapse = " + "))
fm <- paste(fm, "spatial(x_gk, y_gk)", sep = "+")
fm <- as.formula(fm)
sp <- gamboost(fm, data = birds, family = Poisson())

indx <- which(sp$ensemble == 5)
sfit <- attr(sp$data$x[[5]], "fitted")
spat <- rowSums(sapply(indx, function(i) 
    sfit(sp$ensembless[[i]]) * sp$control$nu))
library("akima")
persp(interp(x = birds$x_gk, y = birds$y_gk, z = exp(spat)),
    theta = -40, d = 200, ticktype = "det")
