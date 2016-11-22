#' @include Class-BootPos.R
NULL

################################################################################
#' Class for Nonoverlapping Blocks Bootstrap implementation.
#'
#' \code{NonoverlappingBlocks} is an S4 class that implements the
#' nonoverlapping blocks bootstrap described in Carlstein (1986).
#'
#' \code{NonoverlappingBlocks} extends the S4 class
#' \code{\link{BootPos}} and the remarks made in its documentation
#' apply here as well.
#'
#' The Nonoverlapping Blocks Bootstrap method of Carlstein (1986) resamples blocks
#' randomly, with replacement from the collection of nonoverlapping blocks of
#' length \code{l} that start with observation \eqn{(i-1)l + 1}, where 
#' \eqn{i = 1, \ldots \lfloor N / l \rfloor}.
#' A more precise description of the procedure can also be found in
#' Lahiri (1999), p. 389.
#'
#' @name   NonoverlappingBlocks-class
#' @aliases NonoverlappingBlocks
#' @exportClass NonoverlappingBlocks
#'
#' @keywords S4-classes
#'
#' @seealso \code{\link{getPositions-NonoverlappingBlocks}}
#'
#' @references
#' Carlstein, E. (1986). The use of subseries methods for estimating the
#' variance of a general statistic from a stationary time series.
#' \emph{The Annals of Statistics}, \bold{14}, 1171--1179.
################################################################################

setClass(
    Class = "NonoverlappingBlocks",
    contains = "BootPos"
)

setMethod(
    f = "initialize",
    signature = "NonoverlappingBlocks",
    definition = function(.Object, l, N) {

      .Object@l <- l
      .Object@N <- N

      # Return object
      return(.Object)
    }
)

################################################################################
#' Get Positions for the Nonoverlapping Blocks Bootstrap.
#'
#' @name getPositions-NonoverlappingBlocks
#' @aliases getPositions,NonoverlappingBlocks-method
#'
#' @param object a \code{NonoverlappingBlocks} object; used to specify the parameters
#'                \code{N}, \code{l} and the type of the bootstrap.
#' @param B Number of independent repetitions to bootstrap.
#'
#' @return a matrix of dimension \code{[N,B]} where each column gives the
#'         positions in which to reorder the observations to yield one
#'          bootstrap replication.
################################################################################

setMethod(f = "getPositions",
    signature = "NonoverlappingBlocks",
    definition = function(object, B=1) {

    N <- object@N
    l <- object@l
    nBlocks <- ceiling(N/l)

    positions <- c()

    for (b in 1:B) {
      blocks <- matrix(ncol=nBlocks, nrow=l)
      blocks[1,] <- floor(runif(n=nBlocks, min=0,max=floor(N/l)))*l+1
      if (l > 1) {
        for (i in 2:l) {
          blocks[i,] <- blocks[1,]+i-1
        }
      }
      positions <- c(positions,as.vector(blocks)[1:N])
    }

    return(matrix(positions,nrow=N))
  }
)


################################################################################
#' Create an instance of the \code{\link{NonoverlappingBlocks}} class.
#'
#' @name NonoverlappingBlocks-constructor
#' @aliases nonoverlappingBlocks
#' @export
#'
#' @keywords Constructors
#'
#' @param l the block length for the block bootstrap methods
#' @param N number of available observations to bootstrap from
#'
#' @return Returns an instance of \code{NonoverlappingBlocks}.
################################################################################

nonoverlappingBlocks <- function( l, N ) {

  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }

  obj <- new(
      Class = "NonoverlappingBlocks",
      l = l,
      N = N
  )

  return(obj)
}
