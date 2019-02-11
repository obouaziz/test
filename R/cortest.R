#' Calibrated tests for correlation between paired samples
#'
#' Tests the association/correlation for continuous paired samples using corrected versions of the Pearson correlation test and Kendall tau test. These two tests are asymptotically calibrated.
#' @param x,y the two continuous variables. Must be of same length.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' @param method a character string indicating which test to implement. Can be \code{"pearson"} or \code{"kendall"}.
#' @details Two tests are implemented. The null hypothesis for the corrected Pearson test is: H0 Cor(X,Y)=0 where Cor represents the Pearson correlation coefficient.
#' The alternative is specified by the \code{alternative} argument. The null hypothesis for the corrected Kendall test is: H0 tau=0 where rho represents the Kendall's tau coefficient.
#' Both tests are asymptotically calibrated in the sense that the rejection probability under the null hypothesis is asymptotically equal to 5\%.
#' @return Returns the result of the test with its corresponding p-value, the value of the test statistic and the estimated value of the Pearson correlation coefficient
#' or Kendall's tau.
#' @keywords test
#' @export
#' @examples
#' n=100
#' M=10000 #number of replications
#' testone=function(n){
#' X=rnorm(n,0,1)
#' epsi=rnorm(n,0,1)
#' Y=X^2+0.3*epsi
#' list(test1=cortest(X,Y)$p.value,test2=cor.test(X,Y)$p.value) #cor.test is the standard Pearson test
#' }
#' res1=res2=rep(NA,M)
#' # Replications to check if the the corrected Pearson test and the standard test are well calibrated
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05)  #0.0495
#' mean(res2<0.05)  #0.3674
#'
#' #Replications with Kendall test (takes long time to run)
#' M=1000
#' testone=function(n){
#' X=rnorm(n,0,1)
#' epsi=rnorm(n,0,1)
#' Y=X^2+0.3*epsi
#' list(test1=cortest(X,Y)$p.value,test2=cor.test(X,Y)$p.value,test3=cortest(X,Y,method="kendall")$p.value,test4=cor.test(X,Y,method="kendall")$p.value) #cor.test is the standard Pearson or Kendall correlation test
#' }
#' res1=res2=res3=res4=rep(NA,M)
#' # Replications to check if the the tests are well calibrated
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' #res3[i]=result$test3
#' #res4[i]=result$test4
#' }
#' mean(res1<0.05)  #0.039
#' mean(res2<0.05)  #0.346
#' #mean(res3<0.05) #0.044
#' #mean(res4<0.05) #0.154
#'

cortest <- function(x,y,...) {UseMethod("cortest")}
#' @export
cortest.default=function(X,Y,alternative="two.sided",method="pearson")
{
  if (length(X)!=length(Y)) stop("'x' and 'y' must have the same length")
  n <- length(X)
  if (method=="pearson"){
  X_cent<-(X-mean(X))
  Y_cent<-(Y-mean(Y))
  R <- X_cent*Y_cent
  Tn <- sqrt(n)*((n-1)/n)*cov(X,Y)/sqrt(((n-1)/n)*var(R))
  estimate=sum(R)/(sqrt(sum(X_cent^2)*sum(Y_cent^2)))
  if (alternative=="two.sided" | alternative=="t"){
  Pval <- 2*(1-pt(abs(Tn),n-2))}
  if (alternative=="less"| alternative=="l"){
    Pval <- pt(Tn,n-2)}
  if (alternative=="greater"| alternative=="g"){
    Pval <- 1-pt(Tn,n-2)}
  }
  if (method=="kendall"){
  R <- array(0,dim=c(n,n))
  S <- array(0,dim=c(n,n))
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      R[i,j]<-((X[j]-X[i])*(Y[j]-Y[i]))>0
      S[i,j]<-(X[j]>X[i]) & (Y[j]>Y[i])
    }
  }
  H <- apply(S,1,mean)+apply(S,2,mean)
  V <- var(H)
  estimate=2*(sum(R)/(n*(n-1))-0.5)
  Tn <- sqrt(n)*estimate/(4*sqrt(V))
  if (alternative=="two.sided" | alternative=="t"){
  Pval <- 2*(1-pnorm(abs(Tn)))}
  if (alternative=="less"| alternative=="l"){
    Pval <- pnorm(Tn)}
  if (alternative=="greater"| alternative=="g"){
    Pval <- 1-pnorm(Tn)}
  }
  result <- list(statistic=Tn, p.value=Pval, estimate=estimate,alternative=alternative,method=method)
  class(result)<-"test"
  return(result)
}

#' @export
print.test <- function(x, ...)
{
  corval=x$estimate
  names(corval)<-"cor"
  if (x$method=="pearson"){
  cat("\nCorrected Pearson correlation test\n\n")}
  if (x$method=="kendall"){
    cat("\nCorrected Kendall correlation test\n\n")}
  cat(paste("t = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
  if (x$alternative=="two.sided" | x$alternative=="t"){
  cat("alternative hypothesis: true correlation is not equal to 0\n")}
  if (x$alternative=="less" | x$alternative=="l"){
    cat("alternative hypothesis: true correlation is less than 0\n")}
  if (x$alternative=="greater" | x$alternative=="g"){
    cat("alternative hypothesis: true correlation is greater than 0\n")}
  cat("sample estimates:\n")
  print(corval)
}




