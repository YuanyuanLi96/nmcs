
<!-- README.md is generated from README.Rmd. Please edit that file -->
# nmcs

<!-- badges: start -->
<!-- badges: end -->
The goal of nmcs is to obtain a nested model confidence set(NMCS) and the LogP uncertainty measure for a given shrinkage model selection method.

## Installation

You can install the released version of nmcs from Github with:

``` r
devtools::install_github("YuanyuanLi96/nmcs")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(nmcs)
## basic example code
set.seed(0)
n=100
p=10
B=200
X = matrix(rnorm(n*p),nrow=n,ncol=p)
true_b = c(1:3, rep(0,p-3))
predy = X %*% true_b
#Gaussian
Y=Generate_Y(predy,sigmasq = 1,n=n)
alpha=c(0.05,0.1,0.3)
result=NMCS(Y, X, alpha=alpha,B=B)
summary_NMCS(result,alpha=alpha)#NMCS result
#> $hat_M
#> [1] 1 2 3
#> 
#> $MCS.frame
#>     CL  bcp width     lbm     ubm
#> 1 0.95 0.97     1    2, 3 1, 2, 3
#> 2 0.90 0.97     1    2, 3 1, 2, 3
#> 3 0.70 0.87     0 1, 2, 3 1, 2, 3
result$hat_logP#LogP measure
#> [1] -2.040221
#Binomial
Y=Generate_Y(predy, n=n, family = "binomial")
result=NMCS(Y, X, family="binomial",alpha=alpha, B=B)
summary_NMCS(result,alpha=alpha)#NMCS result
#> $hat_M
#> [1] 1 2 3
#> 
#> $MCS.frame
#>     CL   bcp width  lbm              ubm
#> 1 0.95 0.950     6      1, 2, 3, 4, 5, 9
#> 2 0.90 0.920     4            1, 2, 3, 5
#> 3 0.70 0.765     2 2, 3       1, 2, 3, 5
result$hat_logP#LogP measure
#> [1] -0.7236064
#GAM
Xn=X
Xn[,2]=-1/3*X[,1]^3+rnorm(n)
predy_n = Xn %*% true_b
Yn=Generate_Y(predy_n, n=n, family = "gam")
result=NMCS(Yn, Xn, family="gam",alpha=alpha, B=B)
summary_NMCS(result,alpha=alpha)#NMCS result
#> $hat_M
#> [1] 1 2 3 5
#> 
#> $MCS.frame
#>     CL   bcp width lbm                     ubm
#> 1 0.95 0.970     8     1, 2, 3, 5, 6, 7, 8, 10
#> 2 0.90 0.920     7        1, 2, 3, 5, 7, 8, 10
#> 3 0.70 0.785     6            1, 2, 3, 5, 7, 8
result$hat_logP#LogP measure
#> [1] -0.02531781
```
