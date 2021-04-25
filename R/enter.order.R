#' Entering order of predictors in a regularized path
#'
#' This function allows you to obtain the entering order of predictors by a given regularized
#' solution path.
#' @param beta a matrix with each column being the norm of coefficients
#' for a given tuning parameter \code{lambda}. The number
#' of columns is \code{length(lambda)}.
#' @param eps tolerance when extracting non-zero coefficients. Default value is
#' \code{1e-6}.
#' @return A vector containing the indexes of predictors following their entering order.
#' @export
#' @examples
#' x = matrix(rnorm(100 * 10), 100, 10)
#' y = rnorm(100)
#' library(glmnet)
#' fit1 = glmnet(x, y)
#' enter.order(fit1$beta)

enter.order = function(beta,eps=1e-6){
  enterall = factor(apply(beta,1,function(x)which(x>eps)[1]))
  #break tie by their value
  variable.order=NULL
  for (l in levels(enterall)){
    vnam= which(enterall==l)
    variable.order= c(variable.order,vnam[order(beta[vnam,as.numeric(l)],decreasing = T)])
  }
  variable.order=c(variable.order, as.numeric(which(is.na(enterall))))
  return(variable.order)
}
