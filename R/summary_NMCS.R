#library(ggplot2)
#library(knitr)
#' Summary of nested model confidence sets
#'
#' This is a generic function used to produce result summaries of nested model confidence sets.
#' @param nmcs.r a result of \code{NMCS}.
#' @param predictors The indexes of all predictors. Default value is \code{1:p}.
#' @param alpha A sequence of significance levels. Default value is 0.05.
#' @return Returns a list including the predictors indexes of selected model, and a
#' dataframe including all model confidence sets for user-given \code{alpha} levels.
#' \item{hat_M}{indexes of predictors in the selected model.}
#'  \item{MCS.frame}{a dataframe containing the information about model confidence sets.
#'  \describe{
#'  \item{CL}{confidence levels.}
#'   \item{bcp}{the bootstrap coverage probabilities.}
#' \item{width}{width of NMCS, which equals the size difference between lbm and ubm.}
#' \item{lbm}{lower bound models.}
#' \item{ubm}{upper bound models.}
#' }
#' }
#' @export
#' @examples
#' set.seed(0)
#' n=50
#' p=10
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' true_b = c(1:3, rep(0,p-3))
#' Y = X %*% true_b + rnorm(n)
#' alpha=c(0.05,0.1)
#' result=NMCS(Y, X, alpha=alpha)
#' summary_NMCS(result,alpha=alpha)


summary_NMCS=function(nmcs.r, alpha, predictors=NULL){
  if(is.null(predictors))predictors=1:length(nmcs.r$hat_M$var.order)
  hat_M= sort(predictors[nmcs.r$hat_M$var_M])
  nmcs.r.var=data.frame(CL=sapply(alpha,function(x)1-x),
                      bcp=sapply(1:length(alpha),function(x)
                        nmcs.r$mcs[[x]]$BCP),
                      width = sapply(1:length(alpha),function(x)
                        nmcs.r$mcs[[x]]$width),
                      lbm=matrix(lapply(1:length(alpha),function(x)
                        sort(predictors[nmcs.r$mcs[[x]]$LB]))),
                      ubm=matrix(lapply(1:length(alpha),function(x)
                        sort(predictors[nmcs.r$mcs[[x]]$UB]))))
  return(list(hat_M=hat_M, MCS.frame=nmcs.r.var))
}

