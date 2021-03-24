#' Simulate the distribution of the test statistic for the robust independence test of Kolmogorov-Smirnov's type
#'
#' For two independent continuous uniform variables on [0,1] compute the maximal distance
#' between the joint empirical cumulative distribution function and the product of the marginal
#' empirical cumulative distribution functions using Monte-Carlo simulations.
#' @param n the size of the sample.
#' @param N the number of replications in the Monte-Carlo simulation.
#' @details Let (x1,y1), ..., (x_n,y_n) be a bivariate sample of \code{n} independent continuous uniform variables.
#' Its corresponding bivariate e.c.d.f. (empirical cumulative distribution function)
#' Fn is defined as:
#'
#' Fn(t1,t2) = #\code{{xi<=t1,yi<=t2}/n = sum_{i=1}^n Indicator(xi<=t1,yi<=t2)/n}.
#'
#' Let Fn(t1) and Fn(t2) be the marginals e.c.d.f. Based on N Monte_Carlo simulations,
#' the function computes the e.c.d.f. of
#'
#' \code{n^(1/2) sup_{t1,t2} |Fn(t1,t2)-Fn(t1)*Fn(t2)|}.
#'
#' @return Returns the e.c.d.f. based on the N Monte_Carlo simulations. The returned object
#' is a stepfun object obtained from the function \code{ecdf}.
#' @seealso \code{\link{robustest}}, \code{\link{stat_robustest}}, \code{\link{ecdf2D}}.
#' @export

simulecdf<-function(n,N){
  max_val<-rep(NA, length=N)
  for(i in 1:N)
  {
    x<-runif(n)
    y<-runif(n)
    max_val[i]<-stat_robustest(x,y)
  }
  result<-ecdf(max_val)
  return(result)
}

