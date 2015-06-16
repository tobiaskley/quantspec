#' @include Class-BootPos.R
NULL

################################################################################
#' Class for Circular Blocks Bootstrap implementation.
#'
#' \code{CircularBlocks} is an S4 class that implements the moving blocks
#' bootstrap described in Politis and Romano (1992).
#'
#' \code{CircularBlocks} extends the S4 class
#' \code{\link{BootPos}} and the remarks made in its documentation
#' apply here as well.
#'
#' The Circular Blocks Bootstrap method of Politis and Romano (1992) resamples blocks
#' randomly, with replacement from the collection of overlapping blocks of
#' length \code{l} that start with observation 1, 2, \ldots, \code{N}.
#' In contrast to \code{\link{MovingBlocks}} and \code{\link{NonoverlappingBlocks}},
#' the Circular Blocks Bootstrap uses elements from the periodically extended series.
#' A more precise description of the procedure can also be found in
#' Lahiri (1999), p. 389.
#'
#' @name   CircularBlocks-class
#' @aliases CircularBlocks
#' @exportClass CircularBlocks
#'
#' @keywords S4-classes
#'
#' @seealso \code{\link{getPositions-CircularBlocks}}
#'
#' @references
#' Politis, D. and Romano, J. P. (1992). A circular block resampling procedure
#' for stationary data. In \emph{Exploring the Limits of Bootstrap}
#' (R. Lepage and L. Billar, eds.), 263--270. Wiley, New York.
################################################################################

setClass(
    Class = "CircularBlocks",
    contains = "BootPos"
)

setMethod(
    f = "initialize",
    signature = "CircularBlocks",
    definition = function(.Object, l, N) {

      .Object@l <- l
      .Object@N <- N

      # Return object
      return(.Object)
    }
)

################################################################################
#' Get Positions for the Circular Blocks Bootstrap.
#'
#' @name getPositions-CircularBlocks
#' @aliases getPositions,CircularBlocks-method
#'
#' @param object a \code{CircularBlocks} object; used to specify the parameters
#'                \code{N}, \code{l} and the type of the bootstrap.
#' @param B Number of independent repetitions to bootstrap.
#'
#' @return a matrix of dimension \code{[N,B]} where each column gives the
#'         positions in which to reorder the observations to yield one
#'          bootstrap replication.
################################################################################

setMethod(f = "getPositions",
    signature = "CircularBlocks",
    definition = function(object, B=1) {

    N <- object@N
    l <- object@l
    
    res <- matrix(nrow=N, ncol=B)

    for (b in 1:B) {
      r <- 1
      while (r <= N) {
        res[r, b] <- floor(runif(n=1, min=1,max=N+1))
        if (l > 1) {
          for (i in 2:l) {
            if (r+i-1 <= N) {
              res[r+i-1, b] <- (res[r, b]+i-2) %% N + 1
            }
          }
        }
        r <- r + l
      }
    }

    return(res)
  }
)


################################################################################
#' Create an instance of the \code{\link{CircularBlocks}} class.
#'
#' @name CircularBlocks-constructor
#' @aliases circularBlocks
#' @export
#'
#' @keywords Constructors
#'
#' @param l the block length for the block bootstrap methods
#' @param N number of available observations to bootstrap from
#'
#' @return Returns an instance of \code{CircularBlocks}.
################################################################################

circularBlocks <- function( l, N ) {

  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }

  obj <- new(
      Class = "CircularBlocks",
      l = l,
      N = N
  )

  return(obj)
}
