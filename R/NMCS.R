##################################Main function
#' Nested Model Confidence Set and LogP measure
#'
#' This function allows you to obtain a nested model confidence set and the LogP uncertainty measure
#' for a given shrinkage model selection method.
#' @param Y response variable.
#' @param X covariates matrix, of dimension nobs \eqn{\times} nvars;each row is an observation vector.
#' @param family response type. Either a character string representing one of
#' the families: \code{"gaussian"}, \code{"binomial"}, \code{"gam"}, \code{"poission"}.
#' Default is \code{"gaussian"}.
#' @param B number of bootstrap replicates to perform; Default value is 200.
#' @param alpha Significance level(s). The confidence level of NMCS set is 1-\code{alpha}.
#' Default value is 0.05.
#' @param delta A small positive number added inside of LogP when the bootstrap
#' probability of selected model is 1.
#' @param penalty Default value is \code{"adlasso"}; user can choose from  \code{"adlasso"},
#'   \code{"lasso"},  \code{"scad"}.
#' @param tune method of tuning parameter \eqn{\lambda}. Default method is
#'  \code{"bic"}; user can choose from  \code{"bic"},
#'  \code{"aic"},  \code{"cv"}(stands for "cross validation").
#' @return The NMCS method returns an object of class “NMCS”. An object of
#' class “NMCS” is a list containing at least the following components:
#' \item{mcs}{a list containing \code{alpha} level, and the Bootstrap coverage probability,
#' width, lower bound model, upper bound model
#' of corresponding \code{(1-alpha)%} confidence set.}
#' \item{hat_prob}{the Bootstrap probability for single selected model.}
#' \item{hat_logp}{the LogP measure.}
#' \item{hat_M}{a list containing all the information about the selected model
#' based on original data.
#' \describe{
#' \item{len}{Size of the selected model, which equals the number of non-zero
#'      coefficients.}
#' \item{var.order}{Entering order of all predictors.}
#' \item{beta}{the estimated coefficients of the selected model.}
#' \item{predy}{the fitted values by a linear predictor \eqn{\eta=X\beta}.}
#' \item{sigmasq}{Mean sum of residual squares.}
#' }
#' }
#' @keywords NMCS LogP
#' @export
#' @examples
#' n=100
#' p=10
#' B=200
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' true_b = c(1:3, rep(0,p-3))
#' predy = X %*% true_b
#' #Gaussian
#' Y=Generate_Y(predy,sigmasq = 1,n=n)
#' alpha=c(0.05,0.1,0.3)
#' result=NMCS(Y, X, alpha=alpha,B=B)
#' output_NMCS(result,alpha=alpha)#NMCS result
#' result$hat_logP#LogP measure
#' #Binomial
#' Y=Generate_Y(predy, n=n, family = "binomial")
#' result=NMCS(Y, X, family="binomial",alpha=alpha, B=B)
#' output_NMCS(result,alpha=alpha)#NMCS result
#' result$hat_logP#LogP measure
#' #GAM
#' Xn=X
#' Xn[,2]=-1/3*X[,1]^3+rnorm(n)
#' predy_n = Xn %*% true_b
#' Yn=Generate_Y(predy_n, n=n, family = "gam")
#' result=NMCS(Yn, Xn, family="gam",alpha=alpha, B=B)
#' output_NMCS(result,alpha=alpha)#NMCS result
#' result$hat_logP#LogP measure



NMCS =function(Y, X, family="gaussian", B=200, alpha=0.05, delta=1e-4,
               penalty="adlasso",tune="bic"){
  n = nrow(X)
  p = ncol(X)#number of predictors exclude intercept

  #select model using original data
  hat_M = sel.method(y=Y, x=X, family=family, penalty,tune)#a list including \hat{M}, \hat{\psi} and rank
  order_M = hat_M$var.order
  sigmasq_M = hat_M$sigmasq
  len_M =hat_M$len
  var_M = hat_M$var_M
  predy_M =hat_M$predy

  ##bootstrap procedure
  order_boot = matrix(0, nrow = B, ncol = p)
  len_boot = rep(0,B)
  for (b in 1:B){
    #generate data under hat_M
    y_b <- Generate_Y(predy_M, sigmasq_M, n, family)

    #select model using bootstrap data
    hat_M_boot = sel.method(y=y_b, x=X,family=family, penalty,tune)#a list including \hat{M}, \hat{\psi} and rank
    order_boot[b,] = hat_M_boot$var.order
    len_boot[b] = hat_M_boot$len
  }

  result =list()
  ##confidence sets
  Conf.I = Conf(var_M, len_boot, order_boot, p, B)
  result$mcs = lapply(alpha, function(x)MCS_func(Conf.I, x,order_M, len_M,p))
  ##prob of \hat{M}_b = \hat{M}
  result$hat_prob =  Conf.I[1,1]
  result$hat_logP = logp_delta(result$hat_prob, delta)
  result$hat_M= hat_M
  return(result)
}




