#' Confidence interval for the median in the one sample and in the paired two sample
#'
#' In the one sample case, compute the confidence interval for the median of the random variable. In the paired two sample case, compute the confidence interval for the median of the
#' difference between the two random variables.
#' @param x,y two continuous variables.
#' @param paired a logical value. If it equals TRUE, you must provide values for \code{x} and \code{y}
#' and the paired test is implemented. If it equals FALSE, only \code{x} must be provided.
#' @param conf.level confidence level for the confidence interval.
#' @details Provide the confidence interval for \code{Med(X)} in the one sample case and for \code{Med(X-Y)} in the two sample case.
#' The confidence interval is based on the rank statistic.
#' @note The paired median confidence interval can be implemented by providing the variables \code{x} and \code{y} or by just providing
#' one vector equal to the difference between \code{x} and \code{y}.
#' @return Returns the confidence interval along with the estimate of the median.
#' @keywords test
#' @seealso mediantest
#' @export
#' @examples
#' #Simulations
#' n=100
#' M=2000 #number of replications
#' res1=res2=res3=rep(NA,M)
#' testone=function(n){
#' D=rchisq(n,df=4)-qchisq(df=4, p=0.5)
#' result=medianCI(D)
#' list(test1=mediantest(D)$p.value,test2=c(result$CI))
#' }
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2[1]
#' res3[i]=result$test2[2]
#' }
#' mean(res1<0.05) #0.049
#' (sum(res2>0)+sum(res3<0))/M #0.0525

medianCI <- function(x,y=NULL,paired=FALSE,conf.level=0.95) {UseMethod("medianCI")}
#' @export
medianCI.default=function(X,Y=NULL,paired=FALSE,conf.level=0.95)
{
  alpha=1-conf.level
  if (paired==TRUE)
  {
    if (is.null(Y)){ stop("'y' is missing for paired sample")}
    if (is.null(X)){ stop("'x' is missing for paired sample")}
    #Perform the paired two sample test
    X <- X-Y
    pairedTest<-TRUE
  }
  if (paired==FALSE)
  {
    if (is.null(X)) stop("'x' is missing")
    if (is.null(Y)==FALSE) stop("there should be no 'y' for the one sample median confidence interval")
    pairedTest<-FALSE
  }
  n <- length(X)
  Xsort=sort(X)
  CIl=Xsort[floor(-qnorm(1-alpha/2)*sqrt(n)/2+n/2)]
  CIr=Xsort[floor(qnorm(1-alpha/2)*sqrt(n)/2+(n+1)/2)]
  result <- list(estimate=median(X),CI=c(CIl,CIr),conf.level=conf.level,pairedTest=pairedTest)
  class(result)<-"medianCI"
  return(result)
}

#' @export
print.medianCI <- function(x, ...)
{
  medval=x$estimate
  names(medval)<-"median"
  if (x$pairedTest==FALSE){
    cat("\nConfidence interval for the median\n\n")
    cat(paste(x$conf.level," % confidence interval:","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
    cat("median estimate:\n")
    print(medval)
  }
  if (x$pairedTest==TRUE){
    cat("\nConfidence interval for the median of the difference\n\n")
    cat(paste(x$conf.level," % confidence interval:","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
    cat("median estimate:\n")
    print(medval)
  }
}




