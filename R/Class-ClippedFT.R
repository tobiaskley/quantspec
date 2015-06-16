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
#' a time series as described in Kley et. al (2015+).
#'
#' For each frequency \eqn{\omega} from \code{frequencies} and level \code{q}
#' from \code{levels} the statistic
#' \deqn{\sum_{t=0}^{n-1} I\{Y_t \leq q\} \mbox{e}^{-\mbox{i} \omega t}}
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
#' @references
#' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2015+).
#' Quantile Spectral Processes: Asymptotic Analysis and Inference.
#' \emph{Bernoulli}, \bold{forthcoming}.
#' [cf. \url{http://arxiv.org/abs/1401.8104}]
#'
#' @seealso
#' For an example see \code{\link{FreqRep}}.
################################################################################
setClass(
    Class = "ClippedFT",
    contains = "FreqRep"
)

setMethod(
    f = "initialize",
    signature = "ClippedFT",
    definition = function(.Object, Y, isRankBased, levels, frequencies, positions.boot, B, resampleEcdf) {

      .Object@Y <- Y
      .Object@isRankBased <- isRankBased
      .Object@levels <- levels
      .Object@frequencies <- frequencies
      .Object@positions.boot <- positions.boot
      .Object@B <- B
      .Object@resampleEcdf <- resampleEcdf

      # Define variables with dimensions
      T <- length(Y)
      K <- length(levels)
      J <- length(frequencies)
      
      levels <- sort(levels)

      # Define a matrix to store I{Y_t <= q_k}
      IndMatrix <- matrix(0, nrow=T, ncol=K*(B+1))
      
      if (B > 0) {
        #pos.boot <- array(getPositions(positions.boot,B), dim = c(T,B))
        pos.boot <- getPositions(positions.boot, B)
      }
      
      # for all b -- ==0 non-bootstrap, >0 bootstrap
      for (b in 0:B) {
        
        if (b > 0 && resampleEcdf) {
          pos <- pos.boot[,b]
        } else {
          pos <- 1:T
        }
        
        # Convert Y to "pseudo data", if isRankBased == TRUE
        if (isRankBased) {
          data <- rank(Y[pos], ties.method = "max") / T
        } else {
          data <- Y[pos]
        }
  
        sortedData <- sort(data)
  
        if ( (b == 0) || (b > 0 && resampleEcdf) ) {
          # Fill the matrix
          t <- 1
          for (i in 1:K) {
            while (t <= T && sortedData[t] <= levels[i]) {t <- t+1}
            if (t > 1) {
              IndMatrix[1:(t-1), b*K + i] <- 1
            }
          }
          IndMatrix[, b*K + 1:K] <- IndMatrix[rank(data), b*K + 1:K]
        } else if ( b > 0 && !resampleEcdf ) {
          IndMatrix[, (b*K+1):((b+1)*K)] <- IndMatrix[pos.boot[,b], 1:K]
        }
      }
      cfft <- mvfft(IndMatrix)

      # Modify object to return (only requested frequencies!)
      .Object@values <- array(cfft[unique(T/(2*pi)*frequencies)+1,], dim=c(J,K,B+1))

      # Return object
      return(.Object)
    }
)

################################################################################
#' Create an instance of the \code{\link{ClippedFT}} class.
#'
#' The parameter \code{type.boot} can be set to choose a block bootstrapping
#' procedure. If \code{"none"} is chosen, a moving blocks bootstrap with
#' \code{l=length(Y)} and \code{N=length(Y)} would be done. Note that in that
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
#' @param Y A \code{vector} of real numbers containing the time series from
#'          which to determine the quantile periodogram or a \code{ts} object
#'          or a \code{zoo} object.
#' @param frequencies A vector containing frequencies at which to determine the
#'                    quantile periodogram.
#' @param levels A vector of length \code{K} containing the levels at which the
#'               \code{\link{ClippedFT}} frequency representation is to be
#'               determined.
#' @param isRankBased If true the time series is first transformed to pseudo
#'                    data [cf. \code{\link{FreqRep}}].
#' @param B number of bootstrap replications
#' @param l (expected) length of blocks
#' @param type.boot A flag to choose a method for the block bootstrap; currently
#'                  two options are implemented: \code{"none"}, \code{"mbb"}
#'                  which means to do a moving blocks bootstrap with \code{B}
#'                  and \code{l} as specified. Further options are \code{"nbb"},
#' 								  which means nonoverlapping blocks bootstrap, \code{"cbb"} which
#' 									means circular bootstrap, and \code{"sb"} which stands for
#' 								  stationary bootstrap. 
#' @param resampleEcdf A flag that indicates whether the ecdf used to compute the pseudo
#' 									  data (if \code{isRankBased==TRUE}) is also determined from
#' 										the block bootstraped observations.
#'
#' @return Returns an instance of \code{ClippedFT}.
#'
#' @seealso
#' For an example see \code{\link{FreqRep}}.
################################################################################
clippedFT <- function( Y,
    frequencies=2*pi/length(Y) * 0:(length(Y)-1),
    levels = 0.5,
    isRankBased=TRUE,
    B = 0,
    l = 0,
    type.boot = c("none","mbb","nbb","cbb","sb"),
    resampleEcdf = FALSE) {

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
  frequencies <- frequenciesValidator(frequencies, length(Y))

  type.boot <- match.arg(type.boot, c("none","mbb","nbb","cbb","sb"))[1]
  switch(type.boot,
      "none" = {
        bootPos <- movingBlocks(length(Y),length(Y))},
      "mbb" = {
        bootPos <- movingBlocks(l,length(Y))},
      "nbb" = {
        bootPos <- nonoverlappingBlocks(l,length(Y))},
      "cbb" = {
        bootPos <- circularBlocks(l,length(Y))},
      "sb" = {
        bootPos <- stationaryBlocks(l,length(Y))}
  )

  freqRep <- new(
      Class = "ClippedFT",
      Y = Y,
      isRankBased = isRankBased,
      levels = sort(levels),
      B = B,
      positions.boot = bootPos,
      frequencies = frequencies,
      resampleEcdf
  )

  return(freqRep)
}
