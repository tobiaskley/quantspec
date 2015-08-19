#' @include generics.R
#' @include Class-LagOperator.R

NULL

################################################################################
#' Class to calculate copula covariances from a time series with given levels. 
#' Calculates for each combination of levels \eqn{(\tau_1,\tau_2)}{(tau1,tau2)} 
#' and for all \eqn{k< \code{maxLag}}{k<maxLag} the copula covariances
#' \eqn{Cov(1_{X_0 < \tau_1},1_{X_k < \tau_2})}{Cov(Ind{X0<tau1},Ind{Xk<tau2})}
#' and writes it to \code{values[k]} from its superclass \code{\link{LagOperator}}.
#' 
#' @name ClippedCov-class
#' @aliases ClippedCov
#'
#' @keywords S4-classes
#'
#'
#'
################################################################################

setClass(
  Class = "ClippedCov",
  contains = "LagOperator"
)

#' @importFrom stats quantile
#' @importFrom stats acf
setMethod( 
  f = "initialize",
  signature = "ClippedCov",
  definition = function(.Object,Y,maxLag,levels.1,levels.2,isRankBased) {
    .Object@maxLag = maxLag
    .Object@levels.1 = levels.1
    .Object@levels.2 = levels.2
    .Object@Y = Y
    
    n = length(Y)
    ln.1 = length(levels.1)
    ln.2 = length(levels.2)
    levels.all = union(levels.1,levels.2)
    ln = length(levels.all)
    
    if(isRankBased){
      if(((max(levels.all) > 1) | (min(levels.all)<0)))
      {stop("all levels must be in [0,1] for a Ranked based estimation")}
      Q = quantile(Y,probs = levels.all)
    }
    
    
    Clipped <- matrix(0, nrow=n, ncol=ln)
    for (l in 1:ln) {
      Clipped[,l] <- (Y <= Q[l]) - levels.all[l] 
    }
    pos.1 = match(levels.1,levels.all)
    pos.2 = match(levels.2,levels.all)
    
    .Object@values = array(acf(Clipped,type="covariance", lag.max = maxLag, plot = FALSE, demean = FALSE)$acf[,pos.1,pos.2],dim = c(n,max(ln.1,1),max(ln.2,1)))
    return(.Object) 
  })

################################################################################
#' Create an instance of the \code{\link{ClippedCov}} class.
#'
#' @name ClippedCov-constructor
#' @aliases clippedCov
#' @export
#'
#' @keywords Constructors
#'
#' @param Y Time series to calculate the copula covariance from
#' @param maxLag maximum lag between observations that should be used
#' @param levels.1 a vector of numerics that determines the level of clipping
#' @param levels.2 a vector of numerics that determines the level of clipping
#' @param isRankBased If true the time series is first transformed to pseudo data.
#' 
#' @return Returns an instance of \code{ClippedCov}.
#'
#' @seealso \code{\link{LagOperator}}
#'
#' @examples
#' ccf <- clippedCov(rnorm(100), maxLag = 10, levels.1 =c(0.1,0.5,0.9))
#' dim(getValues(ccf))
#' #print values for levels (.5,.5)
#' print(getValues(ccf)[,2,2])

################################################################################
  clippedCov <- function(Y,maxLag = length(Y) - 1,levels.1 = c(.5),levels.2 = levels.1,isRankBased = TRUE){
    
  if(!(maxLag < length(Y)))
  {maxLag = length(Y) - 1
   warning("maxLag must be smaller then length of dataset, set to maximum")}
  
  if (!((is.vector(levels.1) && is.numeric(levels.2))&&(is.vector(levels.1) && is.numeric(levels.2)))) {
    stop("'levels' needs to be specified as a vector of real numbers")
  }
  

  obj = new("ClippedCov",Y=Y ,maxLag=maxLag ,levels.1=levels.1 ,levels.2=levels.2,isRankBased = isRankBased)
  return(obj)
}