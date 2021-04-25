#Nested MAC Algorithms
##X do not have intercept
#' @import glmnet
#' @import ncvreg
#' @import gamsel

loglik_lam = function(X, y, beta, family) {
  K = dim(beta)[2]
  n =length(y)
  link = cbind(1, X) %*% beta
  yrep = repmat(y, 1, K)
  if (family == "gaussian")
    return(n*log(apply((yrep - link)^2, 2, mean))/2)
  if (family == "poisson")
    return(apply(exp(link) - yrep * link, 2, sum))
  if (family == "binomial")
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
}


getdf = function(coef.beta,eps=1e-6) {
  df=apply(abs(coef.beta) > eps, 2, sum)+1
  return(df)
}

repmat = function(X, m, n) {
  ## R equivalent of repmat (matlab)
  X = as.matrix(X)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}

#Order of predictors by lasso solution path
#beta is the norm of coefficients matrix
enter.list = function(beta,eps=1e-6){
  enterall = factor(apply(beta,1,function(x)which(x>eps)[1]))
  #break tie by their value
  variable.order=NULL
  for (l in levels(enterall)){
    vnam= which(enterall==l)
    variable.order= c(variable.order,vnam[order(beta[vnam,as.numeric(l)],decreasing = T)])
  }
  variable.order=c(variable.order, which(is.na(enterall)))
  return(variable.order)
}


####################Model selection using different tuning methods
cf_norm=function(dat, nbasis){
  m=nrow(dat)/nbasis
  norm=NULL
  for (i in 1:m){
    data_i=dat[nbasis*(i-1)+1:nbasis,]
    norm=rbind(norm, apply(data_i,2, function(x)sqrt(sum(x^2))))
  }
  return(norm)
}
cv.gam<- function(y, x,family,eps=1e-6){
  n=length(y)
  bases=pseudo.bases(x,degree=6,df=4)
  gamsel.cv=cv.gamsel(x,y,bases=bases)
  alpha_norm = sqrt((gamsel.cv$gamsel.fit$alphas)^2)
  #betas_norm = cf_norm(gamsel.cv$gamsel.fit$betas,6)
  #coef_norm = alpha_norm+betas_norm
  index=gamsel.cv$index.min
  hat_M=which(alpha_norm[,index]>eps)
  predy=predict(gamsel.cv$gamsel.fit,x,index=index,type="response")
  sigmasq_est <- mean((predy-y)^2)
  var.order= enter.list(alpha_norm)
  return(list(len= length(hat_M), var.order=var.order,var_M=hat_M, predy=predy
              ,sigmasq=sigmasq_est))
}


sel.method=function(y, x, family, penalty,tune,eps=1e-6){
  n= length(y)
  if (family=="gam"){
    return(cv.gam(y, x,family))
  }
  if(tune=="cv"){
    if(penalty=="scad"){
      cv.fit = cv.ncvreg(x, y, family = family, penalty = "SCAD",nfolds=5)
      betas= cv.fit$fit$beta[-1,]
      cv.1se.ind = min(which(cv.fit$cve<cv.fit$cve[ cv.fit$min]+cv.fit$cvse[ cv.fit$min]))
      coef.beta = cv.fit$fit$beta[, cv.1se.ind]  # extract coefficients at a single value of lambda, including the intercept
    }else{
      lasso_init <- cv.glmnet(x,y, family=family,nfolds = 10)
      betas= lasso_init$glmnet$beta
      coef.beta =as.numeric(coef(lasso_init,s=lasso_init$lambda.min))
      if (penalty=="adlasso"){
        penalty.factor <- abs(coef.beta[-1]+1/sqrt(n))^(-1)
        adalasso <- cv.glmnet(x,y, family=family, penalty.factor=penalty.factor,nfolds = 10)
        betas= adalasso$glmnet$beta
        coef.beta =as.numeric(coef(adalasso,s=adalasso$lambda.min))
      }
    }
  }else{
    if (tune=="aic")k=2
    if (tune=="bic")k=log(n)
    if(penalty=="scad"){
      reg.fit = ncvreg(x, y, family=family,penalty = "SCAD")
      betas = reg.fit$beta#include intercept
    }else{
      lasso_init <- glmnet(x,y, family=family,nfolds = 10)
      betas= coef(lasso_init)
      #df = lasso_init$df
      if (penalty=="adlasso"){
        dev = loglik_lam(x, y, betas, family=family)
        reg.df = getdf(betas[-1, , drop = FALSE])
        obj = 2*dev + k * reg.df
        lambda.ind = which.min(obj)
        coef.beta = betas[, lambda.ind]
        penalty.factor <- abs(coef.beta[-1]+1/sqrt(n))^(-1)
        adalasso <- glmnet(x,y, family=family, penalty.factor=penalty.factor,nfolds = 10)
        betas= rbind(adalasso$a0,adalasso$beta)
      }
    }
    dev = loglik_lam(x, y, betas, family=family)
    reg.df = getdf(betas[-1, , drop = FALSE])
    obj = 2*dev + k * reg.df
    lambda.ind = which.min(obj)
    coef.beta = betas[, lambda.ind]
    betas=betas[-1,]
  }
  #fit glm
  var.order=enter.list(abs(betas))
  hat_M = which(abs(coef.beta[-1])> eps)
  len = length(hat_M)
  if (len==0){
    fit=glm(y~1, family=family)
  }else{
    Xi = x[,hat_M]
    fit=glm(y~Xi, family=family)
  }
  sigmasq_est = sum(fit$residuals^2)/n
  coef.final=rep(0,dim(x)[2])
  coef.final[hat_M]=fit$coefficients[-1]
  return(list(len=len, var.order=var.order,var_M=hat_M, beta=coef.final,
              predy=predict(fit,type = "link"), sigmasq=sigmasq_est))

}
#prob is an element
logp_delta=function(prob, delta){
  logP=rep(0,length(prob))
  for (i in 1:length(prob)){
    if (isTRUE(all.equal(prob[i],0))){
      logP[i]=log(1-delta)
    }else if (isTRUE(all.equal(prob[i], 1))){
      logP[i]=log(delta)
    }else{
      logP[i]=log(1-prob[i])
    }
  }
  return(logP)
}


##################Algorithms for constructing sequence of model confidence sets
#Calculate the Model confidences bounds and their corresponding coverage rate

trunc_bounds = function(x,p){
  if(x<0)x=0
  if(x>p)x=p
  return(x)
}

#alpha ia a given value
MCS_func <- function(bounds, alpha,order_M, len_M,p,eps=1e-6){
  k = which(bounds[1,]- 1+ alpha > -eps)[1]
  BCP=bounds[1,k]
  w=k-1
  j=bounds[2,k]
  #cat("coverage_prob",result$coverage_prob,"\n")
  lower.loc = trunc_bounds(len_M-w+j,p)
  upper.loc = trunc_bounds(len_M+j,p)
  width = upper.loc- lower.loc
  if(upper.loc==0){
    LB=0
    UB=0
  }else{
    if(lower.loc ==0){
      LB = 0
    }else{
      LB =order_M[1:lower.loc]
    }
    UB = order_M[1:upper.loc]
  }
  return(list(alpha= alpha, BCP=BCP, width = width,LB= LB,UB= UB))
}


coverp=function(var_M, len, order, p, B){
  bool=order %in%var_M
  l=ifelse(!all(bool==T),which(bool==F)[1]-1,p)
  ls=ifelse(!all(bool==F),tail(which(bool==T),1), 0)
  d=max(0,len-l)
  ds=max(0, ls-len)
  dr=d+ds
  CP= lapply(0:(2*p),function(x)rep(0,x+1))
  for (w in dr:(2*p)){
    CP[[w+1]][(ds+1):(w-d+1)]=1/B
  }
  return(CP)
}


#Complexity: O(B)
Conf <- function(var_M, len, order.matrix, p, B){
  CP=coverp(var_M, len[1],order.matrix[1,],p,B)
  for (i in 2:B){
    CP=Map("+",CP,coverp(var_M, len[i],order.matrix[i,],p,B))
  }
  maxCP=sapply(CP, function(x)c(max(x), which.max(x)-1))
  #return 2*2p matrix, rows: CP, j*
  return(maxCP)
}



