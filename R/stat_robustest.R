#' Compute the test statistic for the robust independence test of Kolmogorov-Smirnov's type
#'
#' For two continuous variables compute the maximal distance
#' between the joint empirical cumulative distribution function and the product of the marginal
#' empirical cumulative distribution functions.
#' @param x,y the two continuous variables. Must be of same length.
#' @details Let (x1,y1), ..., (x_n,y_n) be a bivariate sample of \code{n} continuous variables.
#' Its corresponding bivariate e.c.d.f. (empirical cumulative distribution function)
#' Fn is defined as:
#'
#' Fn(t1,t2) = #\code{{xi<=t1,yi<=t2}/n = sum_{i=1}^n Indicator(xi<=t1,yi<=t2)/n}.
#'
#' Let Fn(t1) and Fn(t2) be the marginals e.c.d.f. The function returns the value of:
#'
#' \code{n^(1/2) sup_{t1,t2} |Fn(t1,t2)-Fn(t1)*Fn(t2)|}.
#' @return Returns the test statistic of the robust independent test.
#' @seealso \code{\link{robustest}}, \code{\link{simulecdf}}, \code{\link{ecdf2D}}.
#' @export
#' @examples
#' #Simulated data 1
#' x<-c(0.2, 0.3, 0.1, 0.4)
#' y<-c(0.5, 0.4, 0.05, 0.2)
#' stat_robustest(x,y)
#'
#' #Simulated data 2
#' n<-40
#' x<-rnorm(n)
#' y<-x^2+0.3*rnorm(n)
#' plot(x,y)
#' stat_robustest(x,y)
#'
#' #Application on the Evans dataset
#' data(Evans)
#' with(Evans,stat_robustest(CHL[CDH==1],DBP[CDH==1]))

stat_robustest<-function(x,y){
  if (length(x)!=length(x)) stop("'x' and 'y' must have the same length")
  n <- length(x)
  sort_x <- order(x)
  sort_y <- order(y)
  max_val=0
  for(i in 1:(n-1))
  {
    iter=0
    for(j in 1:(n-1))
    {
      if(y[sort_x[j]] <= y[sort_y[i]]){iter<-iter+1}
      val <-abs(iter/n - (i/n)*(j/n))
      if(val>max_val){max_val<-val}
    }
  }
  result<-sqrt(n)*max_val
  return(result)
}


