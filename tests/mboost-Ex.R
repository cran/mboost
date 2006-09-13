### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign("nameEx", 
       local({
	   s <- "__{must remake R-ex/*.R}__"
           function(new) {
               if(!missing(new)) s <<- new else s
           }
       }),
       pos = "CheckExEnv")
## Add some hooks to label plot pages for base and grid graphics
assign("base_plot_hook",
       function() {
           pp <- par(c("mfg","mfcol","oma","mar"))
           if(all(pp$mfg[1:2] == c(1, pp$mfcol[2]))) {
               outer <- (oma4 <- pp$oma[4]) > 0; mar4 <- pp$mar[4]
               mtext(sprintf("help(\"%s\")", nameEx()), side = 4,
                     line = if(outer)max(1, oma4 - 1) else min(1, mar4 - 1),
              outer = outer, adj = 1, cex = .8, col = "orchid", las=3)
           }
       },
       pos = "CheckExEnv")
assign("grid_plot_hook",
       function() {
           pushViewport(viewport(width=unit(1, "npc") - unit(1, "lines"),
                                 x=0, just="left"))
           grid.text(sprintf("help(\"%s\")", nameEx()),
                     x=unit(1, "npc") + unit(0.5, "lines"),
                     y=unit(0.8, "npc"), rot=90,
                     gp=gpar(col="orchid"))
       },
       pos = "CheckExEnv")
setHook("plot.new",     get("base_plot_hook", pos = "CheckExEnv"))
setHook("persp",        get("base_plot_hook", pos = "CheckExEnv"))
setHook("grid.newpage", get("grid_plot_hook", pos = "CheckExEnv"))
assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   .CheckExEnv <- as.environment("CheckExEnv")
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
       pos = "CheckExEnv")
assign("ptime", proc.time(), pos = "CheckExEnv")
grDevices::postscript("mboost-Ex.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), pos = "CheckExEnv")
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
options(warn = 1)    
library('mboost')

assign(".oldSearch", search(), pos = 'CheckExEnv')
assign(".oldNS", loadedNamespaces(), pos = 'CheckExEnv')
cleanEx(); nameEx("FP");
### * FP

flush(stderr()); flush(stdout())

### Name: FP
### Title: Fractional Polynomials
### Aliases: FP
### Keywords: datagen

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




cleanEx(); nameEx("Family");
### * Family

flush(stderr()); flush(stdout())

### Name: Family
### Title: Gradient Boosting Families
### Aliases: Family AdaExp Binomial GaussClass GaussReg Huber Laplace
###   Poisson CoxPH
### Keywords: models

### ** Examples


    Laplace()

    Family(ngradient = function(y, f) y - f, 
           loss = function(y, f) (y - f)^2,
           name = "My Gauss Variant")




cleanEx(); nameEx("blackboost");
### * blackboost

flush(stderr()); flush(stdout())

### Name: blackboost
### Title: Gradient Boosting with Regression Trees
### Aliases: blackboost blackboost_fit
### Keywords: models regression

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- blackboost(dist ~ speed, data = cars,
                          control = boost_control(mstop = 50))
    cars.gb

    ### plot fit
    plot(dist ~ speed, data = cars)
    lines(cars$speed, predict(cars.gb), col = "red")




cleanEx(); nameEx("bodyfat");
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




cleanEx(); nameEx("boost.family-class");
### * boost.family-class

flush(stderr()); flush(stdout())

### Name: boost_family-class
### Title: Class "boost_family": Gradient Boosting Family
### Aliases: boost_family-class show,boost_family-method
### Keywords: classes

### ** Examples


    Laplace()




cleanEx(); nameEx("cvrisk");
### * cvrisk

flush(stderr()); flush(stdout())

### Name: cvrisk
### Title: Cross-Validation
### Aliases: cvrisk
### Keywords: models regression

### ** Examples


  data("bodyfat", package = "mboost")
  tbodyfat <- bodyfat
 
  indep <- names(tbodyfat)[-2]  
  tbodyfat[indep] <- lapply(bodyfat[indep], function(x)
        x <- x - mean(x)
  )
  
  ### fit linear model to data
  model <- glmboost(DEXfat ~ ., data = tbodyfat, 
                    control = boost_control(mstop = 100))

  ### AIC-based selection of number of boosting iterations
  AIC(model)

  ### 10-fold cross-validation
  n <- nrow(tbodyfat)
  k <- 10
  ntest <- floor(n / k)
  cv10f <- matrix(c(rep(c(rep(0, ntest), rep(1, n)), k - 1), 
                    rep(0, n * k - (k - 1) * (n + ntest))), nrow = n)
  cvm <- cvrisk(model, folds = cv10f)
  print(cvm)
  mstop(cvm)
  plot(cvm)

  ### 25 bootstrap iterations
  bs25 <- rmultinom(25, n, rep(1, n)/n)
  cvm <- cvrisk(model, folds = bs25)
  print(cvm)
  mstop(cvm)

  layout(matrix(1:2, ncol = 2))
  plot(cvm)

  ### there seems to be some nonlinearity involved ...
  blackbox <- blackboost(DEXfat ~ ., data = bodyfat)
  cvtree <- cvrisk(blackbox, folds = bs25)
  plot(cvtree)




cleanEx(); nameEx("gamboost");
### * gamboost

flush(stderr()); flush(stdout())

### Name: gamboost
### Title: Gradient Boosting with Componentwise Smoothing Splines
### Aliases: gamboost gamboost.formula gamboost.matrix gamboost_fit
### Keywords: models nonlinear

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

    ### artificial example: sinus transformation
    x <- sort(runif(100)) * 10
    y <- sin(x) + rnorm(length(x), sd = 0.25)
    plot(x, y)
    ### linear model
    lines(x, fitted(lm(y ~ sin(x) - 1)), col = "red")
    ### GAM
    lines(x, fitted(gamboost(y ~ x - 1, 
                    control = boost_control(mstop = 500))), 
          col = "green")




cleanEx(); nameEx("glmboost");
### * glmboost

flush(stderr()); flush(stdout())

### Name: glmboost
### Title: Gradient Boosting with Componentwise Linear Models
### Aliases: glmboost glmboost.formula glmboost.matrix glmboost_fit
### Keywords: models regression

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




cleanEx(); nameEx("methods");
### * methods

flush(stderr()); flush(stdout())

### Name: methods
### Title: Methods for Gradient Boosting Objects
### Aliases: print.glmboost coef.glmboost print.gamboost AIC.gb predict.gb
###   mstop mstop.gbAIC mstop.gb mstop.cvrisk mstop.blackboost fitted.gb
###   logLik.gb
### Keywords: methods

### ** Examples


    ### a simple two-dimensional example: cars data
    cars.gb <- glmboost(dist ~ speed, data = cars, 
                        control = boost_control(mstop = 2000))
    cars.gb

    ### initial number of boosting iterations
    mstop(cars.gb)

    ### AIC criterion
    aic <- AIC(cars.gb, method = "corrected")
    aic

    ### coefficients for optimal number of boosting iterations
    coef(cars.gb[mstop(aic)])
    plot(cars$dist, predict(cars.gb[mstop(aic)]), 
         ylim = range(cars$dist))
    abline(a = 0, b = 1)




cleanEx(); nameEx("wpbc");
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
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
