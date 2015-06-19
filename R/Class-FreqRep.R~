#' @include generics.R
#' @include aux-functions.R
#' @include Class-QSpecQuantity.R
#' @include Class-BootPos.R
NULL

################################################################################
#' Class for Frequency Representation.
#'
#' \code{FreqRep} is an S4 class that encapsulates, for a time
#' series \eqn{(Y_t)_{t=0,\ldots,n-1}}{Y_0,\dots,Y_{n-1}}, the data structures
#' for the storage of a frequency representation. Examples of such frequency
#' representations include
#' \itemize{
#'   \item the Fourier transformation of the clipped time series
#'         \eqn{(\{I\{Y_t \leq q\})}{(I{Y_t <= q})}, or
#'   \item the weighted \eqn{L_1}{L1}-projection of \eqn{(Y_t)} onto an harmonic
#'         basis.
#' }
#' Examples are realized by implementing a sub-class to
#' \code{FreqRep}.
#' Currently, implementations for the two examples mentioned above are available:
#'     \code{\link{ClippedFT}} and
#'     \code{\link{QRegEstimator}}.
#'
#' It is always an option to base the calculations on the pseudo data
#' \eqn{R_{t,n} / n}{R_tn / n} where \eqn{R_{t,n}}{R_tn} denotes the rank of
#' \eqn{Y_t}{Y_t} among \eqn{(Y_t)_{t=0,\ldots,n-1}}{Y_0,\dots,Y_{n-1}}.
#'
#' To allow for a block bootstrapping procedure a number of \code{B} estimates
#' determined from bootstrap replications of the time series which are yield by
#' use of a \code{\link{BootPos}}-object can be stored on initialization.
#'
#' The data in the frequency domain is stored in the array \code{values}, which
#' has dimensions \code{(J,K,B+1)}, where \code{J} is the number of
#' \code{frequencies}, \code{K} is the number of \code{levels} and \code{B} is
#' the number of bootstrap replications requested on intialization.
#' In particular, \code{values[j,k,1]} corresponds to the time series' frequency
#' representation with \code{frequencies[j]} and \code{levels[k]}, while
#' \code{values[j,k,b+1]} is the for the same, but determined from the
#' \code{b}th block bootstrapped replicate of the time series.
#'
#' @name   FreqRep-class
#' @aliases FreqRep
#' @exportClass FreqRep
#'
#' @keywords S4-classes
#'
#' @slot Y The time series of which the frequency representation is to be determined.
#' @slot frequencies The frequencies for which the frequency representation will
#'                   be determined. On initalization
#'                   \code{\link{frequenciesValidator}} is called, so that it
#'                   will always be a vector of reals from \eqn{[0,\pi]}{[0,pi]}.
#'                   Also, only Fourier frequencies of the form
#'                   \eqn{2\pi j / n}{2 pi j / n} with integers \eqn{j} and \eqn{n}
#'                   the \code{length(Y)} are allowed.
#' @slot levels The levels for which the frequency representation will be
#'              determined. If the flag \code{isRankBased} is set to \code{FALSE},
#'              then it can be any vector of reals. If \code{isRankBased} is set
#'              to \code{TRUE}, then it has to be from \eqn{[0,1]}.
#' @slot values The array holding the determined frequency representation. Use a
#'              \code{getValues} method of the relevant subclass to access it.
#' @slot isRankBased A flag that is \code{FALSE} if the determined \code{values}
#'                   are based on the original time series and \code{TRUE} if it
#'                   is based on the pseudo data as described in the Details
#'                   section of this topic.
#' @slot positions.boot An object of type \code{\link{BootPos}},
#'                      that is used to determine the block bootstrapped
#'                      replicates of the time series.
#' @slot B Number of bootstrap replications to perform.
#'
#' @example
#' inst/examples/FreqRep.R
################################################################################

setClass(
    Class = "FreqRep",
    representation=representation(
        Y = "numeric",
        frequencies = "numeric",
        levels = "numeric",
        values = "array",
        isRankBased = "logical",
        positions.boot = "BootPos",
        B = "numeric"
    )
)

################################################################################
#' Get \code{Y} from a \code{\link{FreqRep}} object.
#'
#' @name getY-FreqRep
#' @aliases getY,FreqRep-method
#'
#' @param object \code{FreqRep} of which to get the \code{Y}
#'
#' @return Returns the attribute \code{Y} that's a slot of \code{object}.
################################################################################
setMethod(f = "getY",
    signature = signature("FreqRep"),
    definition = function(object) {
      return(object@Y)
    }
)

################################################################################
#' Get attribute \code{frequencies} from a \code{\link{FreqRep}}.
#'
#' @name getFrequencies-FreqRep
#' @aliases getFrequencies,FreqRep-method
#'
#' @keywords Access-functions
#'
#' @param object \code{FreqRep} from which to get the
#'         \code{frequencies}.
#'
#' @return Returns the \code{frequencies} attribute, as a vector of real numbers.
################################################################################
setMethod(f = "getFrequencies",
    signature = "FreqRep",
    definition = function(object) {
      return(object@frequencies)
    }
)

################################################################################
#' Get attribute \code{levels} from a \code{\link{FreqRep}}.
#'
#' @name getLevels-FreqRep
#' @aliases getLevels,FreqRep-method
#'
#' @keywords Access-functions
#'
#' @param object \code{FreqRep} from which to get the
#'         \code{levels}.
#'
#' @return Returns the \code{levels} attribute, as a vector of real numbers.
################################################################################
setMethod(f = "getLevels",
    signature = "FreqRep",
    definition = function(object) {
      return(object@levels)
    }
)

################################################################################
#' Get values from a frequency representation.
#'
#' For two vectors \code{frequencies} and \code{levels} the values from an
#' \code{object} of type \code{FreqRep} are returned.
#'
#' The two parameters \code{frequencies} and \code{levels} are expected to be
#' vectors of reals; an error is thrown otherwise. If any of the
#' \code{frequencies} or \code{levels} requested is not available from
#' \code{object} a warning is issued, and the values with frequencies and levels
#' closest to the ones requested are returned. Note that the frequencies are
#' transformed to \eqn{[0,\pi]}{[0,pi]} using \code{\link{frequenciesValidator}}
#' when checking if they are available in \code{object}.
#'
#' The returned array of \code{values} is of dimension \code{[J,K,B+1]},
#' where \code{J=length(frequencies)}, \code{K=length(levels)}, and \code{B}
#' denotes the value stored in slot \code{B} of \code{object}. At position
#' \code{(j,k,b)} the returned value is the one corresponding to
#' \code{frequencies[j]} and \code{levels[k]} that are closest to the
#' \code{frequencies} and \code{levels} available in \code{object};
#' \code{\link{closest.pos}} is used to determine what closest to means.
#'
#' @name getValues-FreqRep
#' @aliases getValues,FreqRep-method
#'
#' @keywords Access-functions
#'
#' @param object \code{FreqRep} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels a vector of levels for which to get the values
#'
#' @return Returns data from the array \code{values} that's a slot of
#'          \code{object}.
#'
#' @examples
#' Y        <- rnorm(32)
#' freq     <- 2*pi*c(0:31)/32
#' levels   <- c(0.25,0.5,0.75)
#' cFT      <- clippedFT(Y, freq, levels)
#' V.all    <- getValues(cFT)
#' V.coarse <- getValues(cFT, frequencies = 2*pi*c(0:15)/16, levels = levels)
#' V.fine   <- getValues(cFT, frequencies = 2*pi*c(0:63)/64, levels = levels)
#' V.part   <- getValues(cFT, frequencies = 2*pi*c(0:16)/32, levels = c(0.25))
################################################################################
setMethod("getValues",
    signature(object="FreqRep"),
    function(object,
        frequencies=2*pi*(0:(length(object@Y)-1))/length(object@Y),
        levels=object@levels) {

    # workaround: default values don't seem to work for generic functions?
    if (!hasArg(frequencies)) {
      frequencies <- 2*pi*(0:(length(object@Y)-1))/length(object@Y)
    }
    if (!hasArg(levels)) {
      levels <- object@levels
    }
    # end: workaround

    if (!(is.vector(frequencies)  && is.numeric(frequencies))) {
      stop("'frequencies' needs to be specified as a vector of real numbers")
    }

    if (!(is.vector(levels) && is.numeric(levels))) {
      stop("'levels' needs to be specified as a vector of real numbers")
    }

    if (object@isRankBased && !(prod(levels >= 0) && prod(levels <=1))) {
      stop("'levels' need to be from [0,1] when isRankBased==TRUE")
    }

    ##############################
    ## (Similar) Code also in Class-SmoothedPG!!!
    ## (Similar) Code also in Class-QuantileSD!!!
    ##############################

    # Transform all frequencies to [0,2pi)
    frequencies <- frequencies %% (2*pi)

    # Create an aux vector with all available frequencies
    oF <- object@frequencies
    f <- frequencies

    # returns TRUE if x c y
    subsetequal.approx <- function(x,y) {
      X <- round(x, .Machine$double.exponent-2)
      Y <- round(y, .Machine$double.exponent-2)
      return(setequal(X,intersect(X,Y)))
    }

    C1 <- subsetequal.approx(f[f <= pi], oF)
    C2 <- subsetequal.approx(f[f > pi], 2*pi - oF[which(oF != 0 & oF != pi)])

    if (!(C1 & C2)) {
      warning("Not all 'values' for 'frequencies' requested were available. 'values' for the next available Fourier frequencies are returned.")
    }

    # Select columns
    c.pos <- closest.pos(object@levels,levels)

    if (!subsetequal.approx(levels, object@levels)) {
      warning("Not all 'values' for 'levels' requested were available. 'values' for the next available level are returned.")
    }

    # Select rows
    r1.pos <- closest.pos(oF, f[f <= pi])
    r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])

    J <- length(frequencies)
    K <- length(levels)
    res <- array(dim=c(J, K, object@B+1))

    if (length(r1.pos) > 0) {
      res[which(f <= pi),,] <- object@values[r1.pos,c.pos,]
    }
    if (length(r2.pos) > 0) {
      res[which(f > pi),,] <- Conj(object@values[r2.pos,c.pos,])
    }

    return(res)
  }
)

################################################################################
#' Get \code{isRankBased} from a \code{\link{FreqRep}} object
#'
#' @name getIsRankBased-FreqRep
#' @aliases getIsRankBased,FreqRep-method
#'
#' @keywords Access-functions
#'
#' @param object \code{FreqRep} of which to get the \code{isRankBased}
#'
#' @return Returns the attribute \code{isRankBased} that's a slot of \code{object}.
################################################################################
setMethod(f = "getIsRankBased",
    signature = signature("FreqRep"),
    definition = function(object) {
      return(object@isRankBased)
    }
)

################################################################################
#' Get \code{B} from a \code{\link{FreqRep}} object.
#'
#' @name getB-FreqRep
#' @aliases getB,FreqRep-method
#'
#' @keywords Access-functions
#'
#' @param object \code{FreqRep} of which to get the \code{B}
#' @return Returns the attribute \code{B} that's a slot of \code{object}.
################################################################################
setMethod(f = "getB",
    signature = signature("FreqRep"),
    definition = function(object) {
      return(object@B)
    }
)

################################################################################
#' Get associated \code{\link{BootPos}} from a
#' \code{\link{FreqRep}}.
#'
#' @name getBootPos-FreqRep
#' @aliases getBootPos,FreqRep-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{FreqRep} from which to get the
#'                \code{\link{BootPos}}.
#' @return Returns the \code{\link{BootPos}} object associated.
################################################################################
setMethod(f = "getBootPos",
    signature = "FreqRep",
    definition = function(object) {
      return(object@positions.boot)
    }
)


setMethod(f = "show",
    signature = "FreqRep",
    definition = function(object) {

      values <- getValues(object, frequencies=object@frequencies)
      J <- length(object@frequencies)
      K <- length(object@levels)
      B <- dim(values)[3]

      cat(paste("\n",class(object)," (J=",J,", K=",K,", B+1=",B,")\n", sep=""))

      if (J <= 7) {
        cat("Frequencies: ", round(object@frequencies,4),"\n")
      } else {
        cat("Frequencies: ", round(object@frequencies[1:4],4),"..",round(object@frequencies[(J-2):J],4),"\n")
      }

      if (K <= 10) {
        cat("Levels     : ", round(object@levels,4),"\n")
      } else {
        cat("Levels     : ", round(object@levels[1:5],4),"..",round(object@levels[(K-4):K],4),"\n")
      }

      cat("\nValues:\n")

      resultMatr <- matrix(values[,,1], nrow=J, ncol=K)
      nrowShow <- min(10,nrow(resultMatr))
      ncolShow <- min(4,ncol(resultMatr))

      res <- apply(resultMatr[1:nrowShow, 1:ncolShow, drop=F],c(1,2),
          function(x){complex(real=round(Re(x),3), imaginary=round(Im(x),3))})
      rownames(res) <- round(object@frequencies,3)[1:nrowShow]
      res <- data.frame(res)
      names(res) <- paste("tau=",round(getLevels(object)[1:ncolShow],4), sep="")

      print(res)
    }
)




################################################################################
#' Plot the values of the \code{\link{FreqRep}}.
#'
#' Creates a \code{K} x \code{2} plot depicting a
#' \code{\link{FreqRep}} object.
#' Each of the \code{K} ``lines'' of subplots shows the frequency representation
#' for one value of \eqn{\tau}{tau}. The real and imaginary part are shown on
#' the left and the right, respectively.
#'
#' @name plot-FreqRep
#' @aliases plot,FreqRep,ANY-method
#' @export
#'
#' @param x  The \code{\link{FreqRep}} to plot.
#' @param ratio quotient of width over height of the subplots; use this
#'               parameter to produce landscape or portrait shaped plots.
#' @param frequencies a set of frequencies for which the values are to be
#'                    plotted.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Plots the \code{\link{FreqRep}} for all
#'          \code{frequencies} and \code{levels} specified.
################################################################################

  setMethod(f = "plot",
      signature = signature(x = "FreqRep"),
      definition = function(x,
          ratio = 2,
          frequencies=2*pi*(1:(floor(length(x@Y)/2)))/length(x@Y),
          levels=x@levels) {

    # workaround: default values don't seem to work for generic functions?
    if (!hasArg(frequencies)) {
      frequencies <- 2*pi*(1:(floor(length(x@Y)/2)))/length(x@Y)
    }
    if (!hasArg(levels)) {
      levels <- x@levels
    }
    # end: workaround

tryCatch({

    TT <- length(x@Y)
    K <- length(levels)
    values <- getValues(x,
        frequencies = frequencies, levels=levels)

    def.par <- par(no.readonly = TRUE) # save default, for resetting...

    # TEST
    p <- K
    M1 <- matrix(1:(3*p),ncol=3, byrow=T)
    M <- rbind(c(0,3*p+1,3*p+2),M1,c(0,rep(3*p+3,2)))
    nf <- layout(M, c(lcm(1),rep(ratio,2)), c(lcm(1),rep(1,p),lcm(1)), TRUE)

    # END TEST

    for (i in 1:K) {
      par(mar=c(0,0,0,0))
      plot.new()
      text(0.5,0.5,substitute(paste(tau,"=",k),list(k=round(levels[i],4))), srt=90)

      par(mar=c(2,2,1,1))
      plot(x=frequencies/(2*pi), y=Re(values[,i,1]),
          type="l", xlab="", ylab="")


      plot(x=frequencies/(2*pi), y=Im(values[,i,1]),
          type="l", xlab="", ylab="")
    }


    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,"Real part")
    plot.new()
    text(0.5,0.5,"Imaginary part")

    plot.new()
    text(0.5,0.5,expression(omega/2*pi))

},  error = function(e) e,
    warning = function(w) w,
    finally = {
      par(def.par)  #- reset to default
    })
  }
)
