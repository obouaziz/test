#' Bidimensional Empirical Cumulative Distribution Function
#'
#' Compute an empirical cumulative distribution function for a bivariate continuous distribution.
#' @param x,y the two continuous variables. Must be of same length.
#' @details The bidimensional e.c.d.f. (empirical cumulative distribution function) Fn is a step function with jumps i/n at observation values, where i is
#' the number of tied observations at that value. What about missing values ?
#'
#' For observations (x1,y1), ..., (x_n,y_n), Fn is defined as
#'
#' Fn(t1,t2)=#\code{{xi<=t1,yi<=t2}/n=1/n sum_{i=1}^n Indicator(xi<=t1,yi<=t2)}
#'
#' @return The result is returned as a matrix of dimension (n*n) where the entry (i,j) corresponds to Fn(xi,yj), i=1, ...,n,
#' j=1, ...,n.
#' @keywords test
#' @family test functions
#' @export
#' @examples
#' #Simulated data #1
#' x<-c(0.2, 0.3, 0.1, 0.4)
#' y<-c(0.5, 0.4, 0.05, 0.2)
#' ecdf2D(x,y)
#' #Simulated data #2
#' n<-40
#' x<-rnorm(n)
#' y<-x^2+0.3*rnorm(n)
#' plot(x,y)
#' ecdf2D(x,y)
#' #Real data

ecdf2D <- function(x,y) {UseMethod("ecdf2D")}
#' @export
ecdf2D.default=function(X,Y)
{
  if (length(X)!=length(Y)) stop("'x' and 'y' must have the same length")
  n<-length(X)
  sort_x <- order(x)
  sort_y <- order(y)
  result=matrix(0,n,n)
  for(i in 1:(n))
  {
    cpt=0
    for(j in 1:(n))
    {
      if(y[sort_x[j]] <= y[sort_y[i]]){cpt<-cpt+1}
      result[i,j] <-cpt/n
    }
  }
  #truc=list(result=result)
  #class(truc)<-"ecdf2D"
  #return(truc)
  return(result)
}

#' #' @export
#' print.ecdf2D <- function(x)
#' {
#'   #cat("bidimensional matrix for ecdf of (x,y):\n")
#'   print(x$result)
#' }

