#' Calibrated Wilcoxon rank sum and signed rank tests
#'
#' Compares the distribution between two random variables by testing if one variable tends to take larger (or smaller) values than the other. The test
#' works for independant and paired variables by using corrected versions of the Wilcoxon (or equivalently Mann-Whitney) one and two-sample tests.
#' @param x,y the two continuous variables.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' @param ties.break the method used to break ties in case there are ties in the x or y vectors. Can be \code{"none"} or \code{"random"}.
#' @param paired a logical value. If it equals TRUE, you must provide values for \code{x} and \code{y}
#' and the paired test is implemented. If it equals FALSE, the paired test is implemented when \code{y} is null and
#' only \code{x} is provided and the two sample test (for independent variables) is implemented when both \code{x} and \code{y} are provided.
#' @details The null hypothesis for the corrected Wilcoxon (Mann-Whitney) test is: H0 Med(X-Y)=0 where Med represents the median.
#' The alternative is specified by the \code{alternative} argument: "\code{greater}" means that Med(X-Y)>0 and "\code{less}"
#'  means that Med(X-Y)<0. The null hypothesis for the paired Wilcoxon (Mann-Whitney) test is: H0 Med(D1+D2)=0 where D1 is the difference
#'  between X1 and Y1 taken on the same pair (same with D2 on a different pair). Both tests are asymptotically calibrated in the sense that the rejection probability under the
#'  null hypothesis is asymptotically equal to the level of the test.
#' @note The function can also be called using formulas: type \code{wilcoxtest(x~y,data)} with x the quantitative variable
#' and y a factor variable with two levels. The option \code{ties.break} handles ties in the Wilcoxon test. If \code{ties.break="none"} the ties are ignored, if
#' \code{ties.break="random"} they are randomly broken. For the Wilcoxon rank sum test the ties between the \code{x} and \code{y} are
#' detected and broken (but the ties inside the \code{x} and \code{y} vectors are not changed). For the signed rank test, the ties in the
#' vector \code{x-y} (or in the \code{x} vector in case \code{y=NULL}) are randomly broken.
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic.
#' @keywords test
#' @family test functions
#' @export
#' @examples
#' #Application on the Evans dataset
#' #Description of this dataset is available in the lbreg package
#' with(Evans,wilcox.test(CHL[CDH==0],CHL[CDH==1]))
#' with(Evans,wilcoxtest(CHL[CDH==0],CHL[CDH==1]))
#' wilcoxtest(CHL~CDH,data=Evans) #using formulas
#' wilcoxtest(CHL~CDH,data=Evans,ties.break="random") #same test were ties are randomly broken
#'
#'
#' #For independent samples
#' n=100
#' M=1000 #number of replications
#' testone=function(n){
#' X=runif(n,-0.5,0.5)
#' Y=rnorm(3*n,0,0.04)
#' list(test1=wilcoxtest(X,Y)$p.value,test2=wilcox.test(X,Y)$p.value) #wilcox.test is the standard Wilcoxon test
#' }
#'
#' #Takes a few seconds to run
#' res1=res2=rep(NA,M)
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05) #0.049
#' mean(res2<0.05) #0.174
#'
#' #For paired samples
#' #We use the value of the median of a Gamma distributed variable with shape parameter equal to 1/5
#' #and scale parameter equal to 1. This value is computed from the command qgamma(shape=1/5, scale=1, 0.5)
#' n=100
#' M=1000 #number of replications
#' testone=function(n){
#' D=rgamma(n,shape=1/10,scale=1)-qgamma(shape=1/5, scale=1, 0.5)/2
#' list(test1=wilcoxtest(D)$p.value,test2=wilcox.test(D)$p.value) #wilcox.test is the standard paired Wilcoxon test
#' }
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05) #0.054
#' mean(res2<0.05) #0.074

wilcoxtest <- function(x,y=NULL,...) {UseMethod("wilcoxtest")}#alternative="two.sided",ties.break="none",paired=FALSE
#' @export
wilcoxtest.default=function(X,Y=NULL,alternative="two.sided",ties.break="none",paired=FALSE)
{
  Message=FALSE
  if (paired==TRUE)
  {
    if (is.null(Y)){ stop("'y' is missing for paired test")}
    if (is.null(X)){ stop("'x' is missing for paired test")}
    #Perform the paired two sample test
    X <- X-Y
    n <- length(X)
    if ((length(X)!=length(unique(X)))){
      if (ties.break=="none") {
        warning("The data contains ties!")}
      if (ties.break=="random") {
        Xsort=sort(X,index.return=TRUE)
        if (sum(diff(Xsort$x)==0)>0) {
          index=which(diff(Xsort$x)==0)#which value should be changed in the ordered sample
          X[Xsort$ix[index]]<-X[Xsort$ix[index]]+runif(length(Xsort$ix[index]),-0.00001,0.00001)
        }
        Message=TRUE
      }
    }
    R <- array(0,dim=c(n,n))
    Diag <- vector(mode = "numeric", length = n)
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        R[i,j]<-(X[j]+X[i]>0)
        Diag[i]=R[i,i]
      }
    }
    H <- apply(R,1,mean)
    V <- var(H)
    Tn <- sqrt(n)*((sum(R)-sum(Diag))/(n*(n-1))-0.5)/(2*sqrt(V))
    pairedTest<-TRUE
  }
  if (paired==FALSE)
  {
    if (is.null(X)) stop("'x' is missing")
    if (is.null(Y))
    {
    #Perform the paired two sample test with X being the difference between the variables in the same pair
    n <- length(X)
    if ((length(X)!=length(unique(X)))){
      if (ties.break=="none") {
        warning("The data contains ties!")}
      if (ties.break=="random") {
        Xsort=sort(X,index.return=TRUE)
        if (sum(diff(Xsort$x)==0)>0) {
          index=which(diff(Xsort$x)==0)#which value should be changed in the ordered sample
          X[Xsort$ix[index]]<-X[Xsort$ix[index]]+runif(length(Xsort$ix[index]),-0.00001,0.00001)
        }
      Message=TRUE
      }
    }
    R <- array(0,dim=c(n,n))
    Diag <- vector(mode = "numeric", length = n)
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        R[i,j]<-(X[j]+X[i]>0)
        Diag[i]=R[i,i]
      }
    }
    H <- apply(R,1,mean)
    V <- var(H)
    Tn <- sqrt(n)*((sum(R)-sum(Diag))/(n*(n-1))-0.5)/(2*sqrt(V))
    pairedTest<-TRUE
    } else {
      #Perform the two sample test with X and Y being two independent variables
      n <- length(X)
      m <- length(Y)
      ties=X%in%Y
      if (sum(ties)!=0){
        if (ties.break=="none") {
        warning("The data contains ties between the two vectors!")}
        if (ties.break=="random") {
        X[ties] <- X[ties]+runif(sum(ties),-0.00001,0.00001)
        Message=TRUE
        }
      }
      R <- array(0,dim=c(n,m))
      for(i in 1:n)
      {
        for(j in 1:m)
        {
          R[i,j]<-(Y[j]>X[i])
        }
      }
      H <- apply(R,1,mean)
      G <- apply(R,2,mean)
      V <- var(H)/n + var(G)/m
      Tn <- (mean(R)-0.5)/sqrt(V)
      pairedTest<-FALSE
      }
  }
  if (alternative=="two.sided" | alternative=="t"){
    Pval <- 2*(1-pnorm(abs(Tn)))}
  if (alternative=="less"| alternative=="l"){
    Pval <- pnorm(Tn)}
  if (alternative=="greater"| alternative=="g"){
    Pval <- 1-pnorm(Tn)}
  result <- list(statistic=Tn, p.value=Pval, alternative=alternative,pairedTest=pairedTest,Message=Message)
  class(result)<-"testW"
  return(result)
}

#' @export
wilcoxtest.formula=function(formula,data=list(),...)
{
  mf <- model.frame(formula=formula, data=data)
  response <- attr(attr(mf, "terms"), "response")
  Fact<-factor(mf[[-response]])
  DATA <- setNames(split(mf[[response]], Fact), c("x", "y"))
  result <- wilcoxtest.default(DATA[[1]],DATA[[2]],...)
  return(result)
}

#' @export
print.testW <- function(x, ...)
{
  if (x$pairedTest==FALSE){
    cat("\nCorrected Wilcoxon rank sum test\n\n")
    cat(paste("W = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (X-Y) is not equal to zero\n")}#X and Y tend to take different values
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (X-Y) is negative\n")}# X tends to be smaller than Y
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (X-Y) is positive\n")}}
  if (x$pairedTest==TRUE){
    cat("\nCorrected Wilcoxon signed rank test\n\n")
  cat(paste("W = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
  if (x$alternative=="two.sided" | x$alternative=="t"){
    cat("alternative hypothesis: median (D1+D2) is not equal to zero\n")}#X and Y tend to take different values
  if (x$alternative=="less" | x$alternative=="l"){
    cat("alternative hypothesis: median (D1+D2) is negative\n")}# X tends to be smaller than Y
  if (x$alternative=="greater" | x$alternative=="g"){
    cat("alternative hypothesis: median (D1+D2) is positive\n")}}#X tends to be larger than Y
  if (x$Message==TRUE) {
    cat("\nTies were detected in the dataset and they were randomly broken")
  }
}




