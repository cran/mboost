### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign(".CheckExEnv", as.environment(2), pos = length(search())) # base
## add some hooks to label plot pages for base and grid graphics
setHook("plot.new", ".newplot.hook")
setHook("persp", ".newplot.hook")
setHook("grid.newpage", ".gridplot.hook")

assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   delayedAssign("T", stop("T used instead of TRUE"),
		  assign.env = .CheckExEnv)
	   delayedAssign("F", stop("F used instead of FALSE"),
		  assign.env = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       env = .CheckExEnv)
assign("..nameEx", "__{must remake R-ex/*.R}__", env = .CheckExEnv) # for now
assign("ptime", proc.time(), env = .CheckExEnv)
grDevices::postscript("mboost-Ex.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), env = .CheckExEnv)
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
options(warn = 1)    
library('mboost')

assign(".oldSearch", search(), env = .CheckExEnv)
assign(".oldNS", loadedNamespaces(), env = .CheckExEnv)
cleanEx(); ..nameEx <- "FP"

### * FP

flush(stderr()); flush(stdout())

### Name: FP
### Title: Fractional Polynomials
### Aliases: FP
### Keywords: models

### ** Examples


    data("bodyfat", package = "mboost")
    tbodyfat <- bodyfat
 
    ### map covariates into [1, 2]
    indep <- names(tbodyfat)[-2]
    tbodyfat[indep] <- lapply(bodyfat[indep], function(x) {
        x <- x - min(x)
        x / max(x) + 1
    })
 
    ### generate formula
    fpfm <- as.formula(paste("DEXfat ~ ", paste("FP(", indep, ")", 
                             collapse = "+")))
    fpfm

    ### fit linear model
    bf_fp <- glmboost(fpfm, data = tbodyfat, 
                      control = boost_control(mstop = 3000))

    ### when to stop
    mstop(aic <- AIC(bf_fp))
    plot(aic)

    ### coefficients
    cf <- coef(bf_fp[mstop(aic)])
    length(cf)
    cf[abs(cf) > 0]




cleanEx(); ..nameEx <- "Family"

### * Family

flush(stderr()); flush(stdout())

### Name: Family
### Title: Gradient Boosting Families
### Aliases: Family AdaExp Binomial GaussClass GaussReg Huber Laplace
###   Poisson
### Keywords: misc

### ** Examples


    Laplace()

    Family(ngradient = function(y, f) y - f, 
           loss = function(y, f) (y - f)^2,
           name = "My Gauss Variant")




cleanEx(); ..nameEx <- "blackboost"

### * blackboost

flush(stderr()); flush(stdout())

### Name: blackboost
### Title: Gradient Boosting with Regression Trees
### Aliases: blackboost blackboost_fit
### Keywords: models

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- blackboost(dist ~ speed, data = cars,
                          control = boost_control(mstop = 50))
    cars.gb

    ### plot fit
    plot(dist ~ speed, data = cars)
    lines(cars$speed, predict(cars.gb), col = "red")




cleanEx(); ..nameEx <- "bodyfat"

### * bodyfat

flush(stderr()); flush(stdout())

### Name: bodyfat
### Title: Prediction of Body Fat by Skinfold Thickness, Circumferences,
###   and Bone Breadths
### Aliases: bodyfat
### Keywords: datasets

### ** Examples


    data("bodyfat", package = "mboost")

    ### final model proposed by Garcia et al. (2005)
    fmod <- lm(DEXfat ~ hipcirc + anthro3a + kneebreadth, data = bodyfat)
    coef(fmod)  




cleanEx(); ..nameEx <- "boost.family-class"

### * boost.family-class

flush(stderr()); flush(stdout())

### Name: boost_family-class
### Title: Class "boost_family": Gradient Boosting Family
### Aliases: boost_family-class show,boost_family-method
### Keywords: classes

### ** Examples


    Laplace()




cleanEx(); ..nameEx <- "gamboost"

### * gamboost

flush(stderr()); flush(stdout())

### Name: gamboost
### Title: Gradient Boosting with Componentwise Smoothing Splines
### Aliases: gamboost gamboost.formula gamboost.matrix gamboost_fit
### Keywords: models

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- gamboost(dist ~ speed, data = cars, dfbase = 4, 
                        control = boost_control(mstop = 50))
    cars.gb
    AIC(cars.gb, method = "corrected")

    ### plot fit for mstop = 1, ..., 50
    plot(dist ~ speed, data = cars)    
    tmp <- sapply(1:mstop(AIC(cars.gb)), function(i)
        lines(cars$speed, predict(cars.gb[i]), col = "red"))          
    lines(cars$speed, predict(smooth.spline(cars$speed, cars$dist),
                              cars$speed)$y, col = "green")




cleanEx(); ..nameEx <- "glmboost"

### * glmboost

flush(stderr()); flush(stdout())

### Name: glmboost
### Title: Gradient Boosting with Componentwise Linear Models
### Aliases: glmboost glmboost.formula glmboost.matrix glmboost_fit
### Keywords: models

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- glmboost(dist ~ speed, data = cars, 
                        control = boost_control(mstop = 5000))
    cars.gb

    ### coefficients should coincide
    coef(cars.gb) + c(cars.gb$offset, 0)
    coef(lm(dist ~ speed, data = cars))

    ### plot fit
    plot(dist ~ speed, data = cars)
    lines(cars$speed, predict(cars.gb), col = "red")

    ### alternative loss function: absolute loss
    cars.gbl <- glmboost(dist ~ speed, data = cars, 
                         control = boost_control(mstop = 5000), 
                         family = Laplace())
    cars.gbl

    coef(cars.gbl) + c(cars.gbl$offset, 0)
    lines(cars$speed, predict(cars.gbl), col = "green")

    ### Huber loss with adaptive choice of delta
    cars.gbh <- glmboost(dist ~ speed, data = cars, 
                         control = boost_control(mstop = 5000), 
                         family = Huber())

    lines(cars$speed, predict(cars.gbh), col = "blue")
    legend("topleft", col = c("red", "green", "blue"), lty = 1,
           legend = c("Gaussian", "Laplace", "Huber"), bty = "n")




cleanEx(); ..nameEx <- "methods"

### * methods

flush(stderr()); flush(stdout())

### Name: methods
### Title: Methods for Gradient Boosting Objects
### Aliases: print.glmboost coef.glmboost print.gamboost AIC.gb predict.gb
###   mstop.gbAIC mstop fitted.gb logLik.gb
### Keywords: methods

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- glmboost(dist ~ speed, data = cars, 
                        control = boost_control(mstop = 2000))
    cars.gb

    ### AIC criterion
    aic <- AIC(cars.gb, method = "corrected")
    aic

    ### coefficients for optimal number of boosting iterations
    coef(cars.gb[mstop(aic)])
    plot(cars$dist, predict(cars.gb[mstop(aic)]), 
         ylim = range(cars$dist))
    abline(a = 0, b = 1)




cleanEx(); ..nameEx <- "wpbc"

### * wpbc

flush(stderr()); flush(stdout())

### Name: wpbc
### Title: Wisconsin Prognostic Breast Cancer Data
### Aliases: wpbc
### Keywords: datasets

### ** Examples


    data("wpbc", package = "mboost")

    ### fit logistic regression model with 100 boosting iterations
    coef(glmboost(status ~ ., data = wpbc[,colnames(wpbc) != "time"], 
                  family = Binomial()))




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", env = .CheckExEnv),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
