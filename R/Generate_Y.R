#' Simulate response variable from a glm family
#'
#' This is a generic function that can simulate new observations from a
#'  distribution in exponential family.
#' @param predy Values of a linear predictor that can be written as \eqn{X\beta}.
#' @param sigmasq Variance of errors for a Gaussian linear model. Default value is 1.
#' @param n sample size
#' @param family response type. A character string representing one of
#' the families: \code{"gaussian"}, \code{"binomial"}, \code{"gam"}, \code{"poission"}.
#' Default is \code{"gaussian"}.
#' @return Returns a vector of length \code{n} with elements drawn from
#' a specified family
#' @keywords glm simulation
#' @export
#' @examples
#' set.seed(0)
#' n=50
#' p=10
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' true_b = c(1:3, rep(0,p-3))
#' predy = X %*% true_b
#' #Simulate obs from Gaussian linear model
#' Generate_Y(predy, n)

Generate_Y <- function(predy, sigmasq=1, n, family="gaussian") {
  if (family %in% c("gaussian","gam")){
    Y = predy + rnorm(n,0, sd = sqrt(sigmasq))
  }
  if (family=="binomial"){
    prob=1/(1+exp(-predy))
    Y = rbinom(n=n,size=1, prob=prob)
  }
  if (family == "poisson"){
    mu <- exp(predy)
    Y <- rpois(n=n, lambda=mu)
  }
  return(Y)
}
