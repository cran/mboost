\name{glmboost}
\alias{glmboost}
\alias{glmboost.formula}
\alias{glmboost.matrix}
\alias{glmboost_fit}
\title{ Gradient Boosting with Componentwise Linear Models }
\description{
  Gradient boosting for optimizing arbitrary loss functions where componentwise
  linear models are utilized as base learners.
}
\usage{
\method{glmboost}{formula}(formula, data = list(), weights = NULL, ...)
\method{glmboost}{matrix}(x, y, weights = NULL, ...)
glmboost_fit(object, family = GaussReg(), 
             control = boost_control(), weights = NULL)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. }
  \item{x}{ design matrix. }
  \item{y}{ vector of responses. }
  \item{object}{ an object of class \code{boost_data}, see \code{\link{boost_dpp}}. }
  \item{family}{ an object of class \code{\link{boost_family-class}}, 
                implementing the negative gradient corresponding
                to the loss function to be optimized, by default, squared error loss
                for continuous responses is used.}
  \item{control}{ an object of class \code{\link{boost_control}}.}
  \item{\dots}{ additional arguments passed to callies. }
}
\details{

  A (generalized) linear model is fitted using a boosting algorithm based on componentwise 
  univariate linear models. The fit, i.e., the regression coefficients, can be
  interpreted in the usual way. The methodology is described in 
  Buhlmann and Yu (2003) and Buhlmann (2006).

  The function \code{glmboost_fit} provides access to the fitting
  procedure without data pre-processing, e.g. for cross-validation.


}
\value{
  An object of class \code{glmboost} with \code{\link{print}}, \code{\link{coef}},
  \code{\link{AIC}} and \code{\link{predict}} methods being available.
}
\references{

  Peter Buhlmann and Bin Yu (2003), 
  Boosting with the L2 loss: regression and classification. 
  \emph{Journal of the American Statistical Association}, \bold{98},
  324--339.

  Peter Buhlmann (2006), Boosting for high-dimensional linear models. 
  \emph{The Annals of Statistics}, \bold{34}.

  Peter Buhlmann and Torsten Hothorn (2006),
  Boosting: A statistical perspective. Submitted manuscript.

}
\examples{

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

}
\keyword{models}
\keyword{regression}
