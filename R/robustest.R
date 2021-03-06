#' Robust independence test for two continuous variables of Kolmogorov-Smirnov's type
#'
#' Test the independence between two continuous variables based on the maximum distance
#' between the joint empirical cumulative distribution function and the product of the marginal
#' empirical cumulative distribution functions.
#' @param x,y the two continuous variables. Must be of same length.
#' @param N the number of Monte-Carlo replications if simu=TRUE.
#' @param simu if TRUE a Monte-Carlo simulation with \code{N} replications is used to determine the
#' distribution of the test statistic under the null hypothesis. If FALSE, pre computed tables are used (see Details
#' for more information).
#' @param ties.break the method used to break ties in case there are ties in the x or y vectors. Can be \code{"none"} or \code{"random"}.
#' @details For two continuous variables, robustest tests H0 X and Y are independent
#' against H1 X and Y are not independent.
#'
#' For observations (x1,y1), ..., (x_n,y_n), the bivariate e.c.d.f.
#' (empirical cumulative distribution function) Fn is defined as:
#'
#' Fn(t1,t2) = #\code{{xi<=t1,yi<=t2}/n = sum_{i=1}^n Indicator(xi<=t1,yi<=t2)/n}.
#'
#' Let Fn(t1) and Fn(t2) be the marginals e.c.d.f. The test statistic is defined as:
#'
#' \code{n^(1/2) sup_{t1,t2} |Fn(t1,t2)-Fn(t1)*Fn(t2)|}.
#'
#'Under H0 the distribution of the test statistic is free and is equivalent to
#'the same test statistic computed for two independent continuous uniform variables in [0,1],
#'where the supremum is taken for t1,t2 in [0,1]. Using this result, the distribution of the test
#'statistic is obtained using Monte-Carlo simulations. The user can either use the argument simu=TRUE to
#'perform the Monte-Carlo simulation (with N the number of replications) or simply use the available tables
#'by choosing simu=FALSE. In the latter case, the exact distribution is computed for n=1, ...,150. For 151<=n<=175, the
#'distribution with n=150 is used. For 176<=n<=250, the distribution with n=200 is used.
#'For 251<=n<=400, the distribution with n=300 is used. For 401<=n<=750, the distribution with n=500 is used.
#'For n>=751, the distribution with n=1000 is used. Those tables were computed using 1e^5 replications.
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic.
#' @note Only a two sided alternative is possible with this test. Missing values are removed such that if a value
#' of \code{x} (resp. \code{y}) is missing then the corresponding
#' values of both \code{x} and \code{y} are removed. The test is then implemented on the remaining elements. If \code{ties.break="none"} the ties are ignored, putting
#' mass (nb of ties)/n at tied observations in the computation of the empirical cumulative distribution functions.
#' If \code{ties.break="random"} they are randomly broken.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{vartest}}, \code{\link{mediantest}}, \code{\link{wilcoxtest}}.
#' See also the \code{\link{hoeffd}} function in the \code{\link{Hmisc}} package for the Hoeffding test.
#' @author See \emph{Distribution Free Tests of Independence Based on the Sample Distribution Function}.
#' J. R. Blum, J. Kiefer and M. Rosenblatt, 1961.
#' @export
#' @examples
#' #Simulated data 1
#' x<-c(0.2, 0.3, 0.1, 0.4)
#' y<-c(0.5, 0.4, 0.05, 0.2)
#' robustest(x,y)
#'
#' #Simulated data 2
#' n<-40
#' x<-rnorm(n)
#' y<-x^2+0.3*rnorm(n)
#' plot(x,y)
#' robustest(x,y)
#'
#' #Application on the Evans dataset
#' #Description of this dataset is available in the lbreg package
#' data(Evans)
#' with(Evans,plot(CHL[CDH==1],DBP[CDH==1]))
#' with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1])) #the standard Pearson test
#' with(Evans,cortest(CHL[CDH==1],DBP[CDH==1])) #the robust Pearson test
#' with(Evans,robustest(CHL[CDH==1],DBP[CDH==1])) #the robust independence test
#' #The robust tests give very different pvalues than the standard Pearson test!

#' @export
robustest <- function(x,y,N=50000,simu=FALSE,ties.break="none") {UseMethod("robustest")}
#' @export
robustest.default<-function (X,Y,N=50000,simu=FALSE,ties.break="none"){
  if (length(X)!=length(Y)) stop("'x' and 'y' must have the same length")
  Message=FALSE
  if (sum(is.na(X))!=0)
  {
    na.ind=which(is.na(X))
    X<-X[-na.ind];Y<-Y[-na.ind]
  }
  if (sum(is.na(Y))!=0)
  {
    na.ind=which(is.na(Y))
    X<-X[-na.ind];Y<-Y[-na.ind]
  }
  n <- length(X)
  if (n<3) stop("length of 'x' and 'y' must be greater than 2")
  dupliX=duplicated(X)
  nb_dupliX=sum(dupliX)
  dupliY=duplicated(Y)
  nb_dupliY=sum(dupliY)
  if ((nb_dupliX+nb_dupliY)!=0){
    if (ties.break=="none") {
      warning("The data contains ties!")}
    if (ties.break=="random") {
      Message=TRUE
      if (nb_dupliX!=0){
      X[dupliX]=X[dupliX]+runif(nb_dupliX,-0.00001,0.00001)}
      if (nb_dupliY!=0){
      Y[dupliY]=Y[dupliY]+runif(nb_dupliY,-0.00001,0.00001)}
    }
  }
  if (simu==TRUE){
    ecdf_fun<-simulecdf(n,N)
    Tn<-stat_robustest(X,Y)
    Pval<-1-ecdf_fun(Tn)
  } else {
    #data(ecdf10.Rdata, envir=environment())#paste(ecdf,n,.Rdata,sep="")
    #load(paste("ecdf",n,".Rdata",sep=""))#Tables/
    if (3<=n & n<=150)
    {
      load(system.file(paste("data/ecdf",n,".RData",sep=""),package="test"))
      #data(list=paste("ecdf",n,sep=""))
    } else {
      if (151<=n & n<=175)
      {
        load(system.file("data/ecdf150.RData",package="test"))
        #data(ecdf150)
      } else {
        if (176<=n & n<=250)
        {
          load(system.file("data/ecdf200.RData",package="test"))
          #data(ecdf200)
        } else {
          if (251<=n & n<=400)
          {
            load(system.file("data/ecdf300.RData",package="test"))
            #data(ecdf300)
          } else {
            if (401<=n & n<=750)
            {
              load(system.file("data/ecdf500.RData",package="test"))
              #data(ecdf500)
            } else {
              if (751<=n)
              {
                load(system.file("data/ecdf1000.RData",package="test"))
                #data(ecdf1000)
              }
            }
          }
        }
      }
    }
    Tn<-stat_robustest(X,Y)
    funstep<-stepfun(x1,c(0,y1))
    Pval<-1-funstep(Tn)
  }
  #Pval<-1-ecdf_fun(Tn)
  result <- list(statistic=Tn, p.value=Pval,message=Message)
  class(result)<-"robustest"
  return(result)
}
#' @export
print.robustest <- function(x, ...)
{
  cat("\nRobust independence test for two continuous variables\n\n")
  cat(paste("t = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
  if (x$message==TRUE) {
    cat("\nTies were detected in the dataset and they were randomly broken")
  }
}

