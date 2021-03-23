#' Robust independence test for two continuous variables of Kolmogorov-Smirnov's type
#'
#' Test the independence between two continuous variables based on the maximum distance
#' between the joint empirical cumulative distribution function and the product of the marginal
#' empirical cumulative distribution functions.
#' @param x,y the two continuous variables. Must be of same length.
#' @details For two continuous variables, robustest tests H0 X and Y are independent
#' against H1 X and Y are not independent.
#'
#' For observations (x1,y1), ..., (x_n,y_n), the bivariate e.c.d.f.
#' (empirical cumulative distribution function) Fn is defined as
#'
#' Fn(t1,t2)=#\code{{xi<=t1,yi<=t2}/n= sum_{i=1}^n Indicator(xi<=t1,yi<=t2)/n}
#'
#' Let Fn(t1) and Fn(t2) be the marginals e.c.d.f. The test statistics is defined as:
#'
#' \code{n^(1/2) sup_{t1,t2} |Fn(t1,t2)-Fn(t1)*Fn(t2)|}
#'
#'Under H0 the distribution of the test statistic is free and is equivalent to
#'the same test statistic computed for two continuous uniform variables on [0,1],
#'where the supremum is taken for t1,t2 on [0,1]. Using this result, the distribution of the test
#'statistic is obtained using Monte-Carlo simulations.
#'
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic.
#' @note Only a two sided alternative is possible with this test.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{vartest}}, \code{\link{mediantest}}, \code{\link{wilcoxtest}}.
#' @export
#' @examples
#' #Simulated data 1
#' x<-c(0.2, 0.3, 0.1, 0.4)
#' y<-c(0.5, 0.4, 0.05, 0.2)
#' stat_robustest(x,y)
#' #Simulated data 2
#' n<-40
#' x<-rnorm(n)
#' y<-x^2+0.3*rnorm(n)
#' plot(x,y)
#' stat_robustest(x,y)
#' #Real data

#' @export
robustest <- function(x,y,N=50000,simu=FALSE) {UseMethod("robustest")}
#' @export
robustest.default<-function (X,Y,N=50000,simu=FALSE){
  if (length(X)!=length(Y)) stop("'x' and 'y' must have the same length")
  n <- length(X)
  if (simu==TRUE){
    ecdf_fun<-simulecdf(n,N)
  } else {
    #data(ecdf10.Rdata, envir=environment())#paste(ecdf,n,.Rdata,sep="")
    load(paste("ecdf",n,".Rdata",sep=""))#Tables/
  }
  Tn<-stat_robustest(X,Y)
  Pval<-1-ecdf_fun(Tn)
  result <- list(statistic=Tn, p.value=Pval)
  class(result)<-"robustest"
  return(result)
}
#' @export
print.robustest <- function(x, ...)
{
  cat("\nRobust independence test for two continuous variables\n\n")
  cat(paste("t = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
}

