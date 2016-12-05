#' @include Class-FreqRep.R
#' @include Class-BootPos.R
NULL

################################################################################
#' Class for Fourier transform of the clipped time series.
#'
#' \code{ClippedFT} is an S4 class that implements the necessary
#' calculations to determine the Fourier transform of the clipped time
#' series. As a subclass to \code{\link{FreqRep}} it inherits
#' slots and methods defined there; it servers as a frequency representation of
#' a time series as described in Kley et. al (2016) for univariate time series
#' and in Barunik & Kley (2015) for multivariate time series.
#'
#' For each frequency \eqn{\omega} from \code{frequencies} and level \code{q}
#' from \code{levels} the statistic
#' \deqn{\sum_{t=0}^{n-1} I\{Y_{t,i} \leq q\} \mbox{e}^{-\mbox{i} \omega t}}
#' is determined and stored to the array \code{values}. Internally the methods
#' \code{\link[stats]{mvfft}} and \code{\link[stats]{fft}} are used to achieve
#' good performance.
#'
#' Note that, all remarks made in the documentation of the super-class
#' \code{\link{FreqRep}} apply.
#'
#' @name ClippedFT-class
#' @aliases ClippedFT
#' @exportClass ClippedFT
#'
#' @keywords S4-classes
#' 
#' @slot multipliers.boot multipliers to be used
#'
#' @references
#' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2016).
#' Quantile Spectral Processes: Asymptotic Analysis and Inference.
#' \emph{Bernoulli}, \bold{22}(3), 1770--1807.
#' [cf. \url{http://arxiv.org/abs/1401.8104}]
#' 
#' Barunik, J. & Kley, T. (2015).
#' Quantile Cross-Spectral Measures of Dependence between Economic Variables.
#' [preprint available from the authors]
#'
#' @seealso
#' For an example see \code{\link{FreqRep}}.
################################################################################
setClass(
    Class = "ClippedFT",
    representation=representation(
        multipliers.boot = "BootMultipliers"
    ),
    contains = "FreqRep"
)

#' @importFrom stats mvfft
setMethod(
    f = "initialize",
    signature = "ClippedFT",
    definition = function(.Object, Y, isRankBased, levels, frequencies, positions.boot, multipliers.boot, B) {
      
      .Object@Y <- Y
      .Object@isRankBased <- isRankBased
      .Object@levels <- levels
      .Object@frequencies <- frequencies
      .Object@positions.boot <- positions.boot
      .Object@multipliers.boot <- multipliers.boot
      .Object@B <- B
      
      # Define variables with dimensions
      T <- dim(Y)[1]
      D <- dim(Y)[2]
      K <- length(levels)
      J <- length(frequencies)

      # Define a matrix to store I{Y_t <= q_k}
      IndMatrix <- matrix(0, nrow=T, ncol=K*D*(B+1))
      
      levels <- sort(levels)
      
      ## Create the Indicator matrix
      
      if (class(multipliers.boot) == "NoneMultipliers") {
      
        ## VERSION 1 - works without multipliers!
        ## This used to be an option to choose!!
        resampleEcdf <- TRUE
  
        if (B > 0) {
          pos.boot <- getPositions(.Object@positions.boot, B)
          mult.boot <- getMultipliers(.Object@multipliers.boot, B)
        }
        
        for (b in 0:B) {
          
          ## If bootstrap with resampled ecdf
          if (b > 0 && resampleEcdf) {
            pos <- pos.boot[, b]
          } else {
            pos <- 1:T
          }
          
          # Convert Y to "pseudo data", if isRankBased == TRUE
          if (isRankBased) {
            data <- apply(Y[pos,, drop=F], 2, rank, ties.method = "max") / T
          } else {
            data <- Y[pos,, drop=F]
          }
          
          for (d in 1:D) {
            sortedData <- sort(data[,d])
            
            if ( (b == 0) || (b > 0 && resampleEcdf) ) {
              # Fill the matrix
              t <- 1
              for (i in 1:K) {
                while (t <= T && sortedData[t] <= levels[i]) {t <- t+1}
                if (t > 1) {
                  IndMatrix[1:(t-1), b*(K*D)+(d-1)*K+i] <- 1
                }
              }
              IndMatrix[, b*(K*D)+(d-1)*K+1:K] <- IndMatrix[rank(data[,d]), b*(K*D)+(d-1)*K+1:K]
            }
          }
          
          if ( b > 0  && !resampleEcdf ) {
            IndMatrix[,(b*(K*D)+1):((b+1)*(K*D))] <- IndMatrix[pos.boot[,b],1:(K*D)]
          }
  
        }
        ## END VERSION 1
        } else {

#      ## VERSION 2 - naive!
#
#      if (B > 0) {
#        pos.boot <- getPositions(.Object@positions.boot, B)
#        mult.boot <- getMultipliers(.Object@multipliers.boot, B)
#      }
#
#      for (b in 0:B) {
#        
#          if (b > 0) {
#            pos <- pos.boot[, b]
#            mult <- mult.boot[, b]          
#          } else {
#            pos <- 1:T
#            mult <- rep(1,T)
#          }  
#        
#        for (d in 1:D) {
#          
#          if (isRankBased) {
#            myEdf <- function(x) {
#              myEdf.simple <- function(x) {return( sum( mult * (Y[pos, d] <= x) ) / T )}
#              return( Vectorize(myEdf.simple)(x) )
#            }
#          } else {
#            myEdf <- function(x) {
#              return( x )
#            }
#          }
#          
#          for (i in 1:K) {
#            for (t in 1:T) {
#              IndMatrix[t, b*(K*D)+(d-1)*K+i] <- mult[t] * ( myEdf(Y[pos[t],d]) <= levels[i] )
#            }
#          }
#        }
#      }
#      ## END VERSION 2

        ## VERSION 3 - implemented via Rcpp
  
        pos.boot <- matrix(0:(T-1), ncol=1)
        mult.boot <- matrix(rep(1,T), ncol=1)
        if (B > 0) {
          pos.boot <- cbind(pos.boot, getPositions(.Object@positions.boot, B) - 1)
          mult.boot <- cbind(mult.boot, getMultipliers(.Object@multipliers.boot, B))
        }
  
        IndMatrix <- .generateIndMatrix(Y, pos.boot, mult.boot, levels, isRankBased)
        ## END VERSION 3

      }

      cfft <- mvfft(IndMatrix)
      
      # Modify object to return (only requested frequencies!)
      .Object@values <- array(cfft[unique(T/(2*pi)*frequencies)+1,], dim=c(J,K,D,B+1))
      .Object@values <- aperm(.Object@values, perm=c(1,3,2,4))
      # Return object
      return(.Object)
    }
)

################################################################################
#' Create an instance of the \code{\link{ClippedFT}} class.
#'
#' The parameter \code{type.boot} can be set to choose a block bootstrapping
#' procedure. If \code{"none"} is chosen, a moving blocks bootstrap with
#' \code{l=lenTS(Y)} and \code{N=lenTS(Y)} would be done. Note that in that
#' case one would also chose \code{B=0} which means that \code{getPositions}
#' would never be called. If \code{B>0} then each bootstrap replication would
#' be the undisturbed time series.
#'
#' @name ClippedFT-constructor
#' @aliases clippedFT
#' @export
#'
#' @keywords Constructors
#'
#' @param Y A \code{matrix} of real numbers containing the time series from
#'          which to determine the quantile periodogram as columns, or a
#' 					\code{ts} object or a \code{zoo} object.
#' @param frequencies A vector containing frequencies at which to determine the
#'                    quantile periodogram.
#' @param levels A vector of length \code{K} containing the levels at which the
#'               \code{\link{ClippedFT}} frequency representation is to be
#'               determined.
#' @param isRankBased If true the time series is first transformed to pseudo
#'                    data [cf. \code{\link{FreqRep}}].
#' @param B number of bootstrap replications
#' @param l (expected) length of blocks or the parameter for the multipliers
#' @param type.boot A flag to choose a method for the block bootstrap; currently
#'                  two options are implemented: \code{"none"}, \code{"mbb"}
#'                  which means to do a moving blocks bootstrap with \code{B}
#'                  and \code{l} as specified. Further options are \code{"nbb"},
#' 								  which means nonoverlapping blocks bootstrap, \code{"cbb"} which
#' 									means circular bootstrap, and \code{"sb"} which stands for
#' 								  stationary bootstrap.
#' @param bootMultipliers \code{\link{BootMultipliers}} object; has to be of type
#'                \code{MovingAverageMultipliers} if \code{type.boot=="mult.ma"}, of
#'                type \code{CovarianceMatrixMultipliers} if \code{type.boot=="mult.cov"}
#' 								and of type \code{NoneMultiplier} otherwise (default).
#'
#' @return Returns an instance of \code{ClippedFT}.
#'
#' @seealso
#' For an example see \code{\link{FreqRep}}.
################################################################################
clippedFT <- function( Y,
    frequencies=2*pi/lenTS(Y) * 0:(lenTS(Y)-1),
    levels = 0.5,
    isRankBased=TRUE,
    B = 0,
    l = 0,
    type.boot = c("none","mbb","nbb","cbb","sb","mult.ma","mult.cov"),
    bootMultipliers = noneMultipliers(lenTS(Y))) {
  
  # Verify if all parameters are valid
  Y <- timeSeriesValidator(Y)
  
  if (!(is.vector(frequencies)  && is.numeric(frequencies))) {
    stop("'frequencies' needs to be specified as a vector of real numbers")
  }
  
  if (!(is.vector(levels) && is.numeric(levels))) {
    stop("'levels' needs to be specified as a vector of real numbers")
  }
  
  if (isRankBased && !(prod(levels >= 0) && prod(levels <=1))) {
    stop("'levels' need to be from [0,1] when isRankBased==TRUE")
  }
  
  # Check validity of frequencies
  frequencies <- frequenciesValidator(frequencies, lenTS(Y))
  
  type.boot <- match.arg(type.boot, c("none", "mbb", "nbb", "cbb", "sb", "mult.ma", "mult.cov"))[1]
  switch(type.boot,
      "none" = {
        bootPos <- movingBlocks(lenTS(Y),lenTS(Y))
        if (class(bootMultipliers) != "NoneMultipliers") {
          stop("bootMultipliers has to be of type NoneMultipliers if type.boot=='none'.")
        }},
      "mbb" = {
        bootPos <- movingBlocks(l,lenTS(Y))
        if (class(bootMultipliers) != "NoneMultipliers") {
          stop("bootMultipliers has to be of type NoneMultipliers if type.boot=='mbb'.")
        }},
      "nbb" = {
        bootPos <- nonoverlappingBlocks(l,lenTS(Y))
        if (class(bootMultipliers) != "NoneMultipliers") {
          stop("bootMultipliers has to be of type NoneMultipliers if type.boot=='nbb'.")
        }},
      "cbb" = {
        bootPos <- circularBlocks(l,lenTS(Y))
        if (class(bootMultipliers) != "NoneMultipliers") {
          stop("bootMultipliers has to be of type NoneMultipliers if type.boot=='cbb'.")
        }},
      "sb" = {
        bootPos <- stationaryBlocks(l,lenTS(Y))
        if (class(bootMultipliers) != "NoneMultipliers") {
          stop("bootMultipliers has to be of type NoneMultipliers if type.boot=='sb'.")
        }},
      "mult.ma" = {
        bootPos <- movingBlocks(lenTS(Y),lenTS(Y))
        if (class(bootMultipliers) != "MovingAverageMultipliers") {
          stop("bootMultipliers has to be of type MovingAverageMultipliers if type.boot=='mult.ma'.")
        }},
      "mult.cov" = {
        bootPos <- movingBlocks(lenTS(Y),lenTS(Y))
        if (class(bootMultipliers) != "CovarianceMatrixMultipliers") {
          stop("bootMultipliers has to be of type CovarianceMatrixMultipliers if type.boot=='mult.cov'.")
        }}
  )
  bootMult <- bootMultipliers
  
  freqRep <- new(
      Class = "ClippedFT",
      Y = Y,
      isRankBased = isRankBased,
      levels = sort(levels),
      B = B,
      positions.boot = bootPos,
      multipliers.boot = bootMultipliers,
      frequencies = frequencies
  )
  
  return(freqRep)
}
