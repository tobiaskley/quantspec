#' @include generics.R
#' @include Class-Weight.R
#' @include Class-LagOperator.R
#' 
NULL
################################################################################
#' Class for a lag window type estimator.
#'
#' For a given time series Y a lag window estimator of the Form
#' \deqn{\hat{f}(\omega) = \sum_{|k|< n-1 } K_n(k) \Gamma(Y_0,Y_k) \exp(-i \omega k)} 
#' will be calculated on initalization. The \code{LagKernelWeight} K_n is determined
#' by the slot \code{weight} and the \code{LagOperator} \eqn{\Gamma(Y_0,Y_k)} is defined 
#' by the slot lagOp.
#' 
#'
#' @name LagEstimator-class
#' @aliases LagEstimator
#' @exportClass LagEstimator
#'
#' @keywords S4-classes
#'
#' @slot Y the time series where the lag estimator will be applied one
#' @slot weight a \code{\link{Weight}} object to be used as lag window
#' @slot lagOp a \code{\link{LagOperator}} object that determines which
#'             kind of bivariate structure should be calculated.
#'
################################################################################
setClass(
  Class = "LagEstimator",
  representation=representation(
    Y = "numeric",
    weight = "Weight",
    lagOp = "LagOperator"   
  ),
  contains = "QSpecQuantity"
)

setMethod(
  f = "initialize",
  signature = "LagEstimator",
  definition = function(.Object, Y,frequencies,lagOp,weight){
   
    levels.1 = getLevels(lagOp,1)
    levels.2 = getLevels(lagOp,2)
    
    .Object@Y = Y
    .Object@frequencies = frequencies
    .Object@levels[[1]] = levels.1
    .Object@levels[[2]] = levels.2
    .Object@weight = weight
    .Object@lagOp = lagOp
    
    levels.all = union(levels.1,levels.2)
    ln1 = length(levels.1)
    ln2 = length(levels.2)
    ln = length(levels.all)
    Q <- Y
    Q <- quantile(Y,probs = levels.all)
    
    
    res <- array(dim = c(length(Y),ln1,ln2))
    
    Fhat = getValues(lagOp,levels.1 = levels.1,levels.2 = levels.2)
    n = length(Fhat)/(ln1*ln2)
    Kernel = getValues(weight)
    if(length(Kernel)<n){
      warning("number of 'weights' is too small, filled up with zeroes")
      Kernel = c(Kernel,rep(0,length(Fhat[,1,1])-length(Kernel)))
    }
    if(length(Kernel)>length(Fhat[,1,1])){
      Kernel = Kernel[1:length(Fhat[,1,1])]
    }
    
    for (l1 in 1:ln1) {
      for (l2 in 1:ln2) {  #\tau_2,\tau_1 can be calculated from \tau_1,\tau_2 
        if((l1 > l2)&&!(levels.1[l1]==levels.2[l2])&&is.element(levels.1[l1],levels.2)&&is.element(levels.2[l2],levels.1)){ 
          res[,l1,l2] <- Conj(res[,l2,l1])
          next
        }
        A1 <- Fhat[,l1,l2]
        A1 <- A1# - levels.1[l1] * levels.2[l2]
        A1 <- A1*Kernel
        
        A2 <- Fhat[,l2,l1]
        A2 <- A2# - levels.1[l2] * levels.2[l1]  				 # same as above!
        A2 <- A2*Kernel
        
        A0 <- Fhat[1,l2,l1]
        A0 <- A0# - levels.1[l2] * levels.2[l1]					 # same as above!
        A0 <- A0*Kernel[1]
        print(A1)
        print(A2)
        print(A0)
        res[,l1,l2] <- 1/(2*pi) * (fft(A1) + Conj(fft(A2)) - A0)
      }
    }
    .Object@values = res
    return(.Object)
  })

  ################################################################################
  #' Create an instance of the \code{LagEstimator} class.
  #'
  #' A \code{LagEstimator} object can be created from \code{numeric}, a \code{ts},
  #' or a \code{zoo} object. Also a \code{\link{LagOperator}} and a 
  #' \code{\link{Weight}} object can be used to create different types of 
  #' estimators.
  #'
  #'
  #' @name LagEstimator-constructor
  #' @aliases lagEstimator
  #' @export
  #'
  #' @keywords Constructors
  #'
  #' @param Y a time series (\code{numeric}, \code{ts}, or \code{zoo} object) or a 
  #'          \code{\link{LagOperator}} from which to determine the \code{LagEstimator}
  #' @param frequencies A vector containing (Fourier-)frequencies at which to determine the
  #'                     smoothed periodogram.
  #' @param levels.1 the first vector of levels for which to compute the LagEstimator
  #' @param levels.2 the second vector of levels for which to compute the LagEstimator               
  #' @param weight Object of type \code{\link{Weight}} to be used for smoothing.
  #' @param type if \code{Y} is a time series, this indicates which LagOperator will be used 
  #' @return Returns an instance of \code{LagEstimator}.
  #'
  #' @examples
  #' Y <- rnorm(100)
  #' levels.1 <- c(0.1,0.5,0.9)
  #' weight <- lagKernelWeight(W = WParzen,  bw = 10, K = length(Y))
  #' lagOp <- clippedCov(Y,levels.1 = levels.1)
  #' lagEst <- lagEstimator(lagOp,weight = weight)
################################################################################
lagEstimator <- function(Y,
                         frequencies=2*pi/length(Y) * 0:(length(Y)-1),
                         levels.1 = .5,
                         levels.2 = levels.1,
                         weight = lagKernelWeight(K = length(Y)-1,bw = 25),
                         type = c("clippedCov")
                         ){
  
  #ToDo Check if params are okay
  #     Transform to Fourier frequencies
  if (class(Y) == "numeric") {
    versConstr <- 1
    Y <- Y
  } else if (class(Y) == "ts") {
    versConstr <- 1
    Y <- Y[1:length(Y)]
  } else if (class(Y) == "zoo") {
    versConstr <- 1
    Y <- coredata(Y)
  } else if ( is(Y,"LagOperator")) {
    LagOp = Y
    versConstr <- 2
    if (!hasArg(frequencies)) {
      Z <- Y@Y
      frequencies <- 2*pi/length(Z) * 0:(length(Z)-1)
    }
    
    if (!hasArg(levels.1)) {
      levels.1 <- getLevels(Y,1)
    }
    if (!hasArg(levels.2)) {
      levels.2 <- getLevels(Y,2)
    }
    
  } else {
    stop("object is neither 'numeric', 'ts', 'zoo', nor 'LagOperator'.")
  }
  if(versConstr == 1){
    if (type == "clippedCov"){
      lagOp = clippedCov(Y,levels.1 = levels.1, levels.2 = levels.2)  
    }
    
  }

   return(new(Class = "LagEstimator",Y = Z,frequencies = frequencies,weight = weight, lagOp = lagOp))
  
  
}
  
  
  
  
