
require("mboost")

set.seed(290875)

### for boosting hat matrix checks
fm <- GaussReg()
fm@offset <- function(y, w) 0

### a simple two-dimensional example from `gamboost.Rd'
data("cars")
cars.gb <- gamboost(dist ~ speed, data = cars, df = 4, family = fm,
                    control = boost_control(mstop = 50))
cars.gb
aic <- AIC(cars.gb, method = "corrected")
aic

### plot fit
plot(dist ~ speed, data = cars)
lines(cars$speed, predict(cars.gb[mstop(AIC(cars.gb))]), col = "red")
lines(cars$speed, predict(smooth.spline(cars$speed, cars$dist), cars$speed)$y, 
      col = "green")

#### check boosting hat matrix and subsetting / predict
stopifnot(isTRUE(all.equal(drop(attr(aic, "hatmat") %*% cars$dist),
                           as.vector(predict(cars.gb)))))
stopifnot(isTRUE(all.equal(drop(attr(AIC(cars.gb[25]), "hatmat") %*% cars$dist),
                           as.vector(predict(cars.gb[25])))))
stopifnot(isTRUE(all.equal(drop(attr(AIC(cars.gb[25]), "hatmat") %*% cars$dist),
                           as.vector(fitted(cars.gb[25])))))

### check boosting hat matrix with multiple independent variables
### and weights
data("bodyfat", package = "mboost")
bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4 - 1
indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
bodyfat[indep] <- lapply(bodyfat[indep], function(x) x - mean(x))
bf_gam <- gamboost(bffm, data = bodyfat, control = boost_control(mstop = 10), 
                   weights = runif(nrow(bodyfat)) * 10)
aic <- AIC(bf_gam)

off <- bf_gam$offset
u <- bf_gam$ustart

stopifnot(isTRUE(all.equal(drop(attr(aic, "hatmat") %*% u + off),
                           as.vector(predict(bf_gam)))))
stopifnot(isTRUE(all.equal(drop(attr(aic, "hatmat") %*% u + off),
                           as.vector(fitted(bf_gam)))))
