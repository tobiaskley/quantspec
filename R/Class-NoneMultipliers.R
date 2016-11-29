#' @include Class-BootMultipliers.R
NULL

################################################################################
#' Class for Multipliers that are 1.
#'
#' \code{NoneMultipliers} is an S4 class that implements the not using multipliers
#'
#' \code{NoneMultipliers} extends the S4 class
#' \code{\link{BootMultipliers}} and the remarks made in its documentation
#' apply here as well.
#'
#' TODO: Add more description here!
#'
#' @name   NoneMultipliers-class
#' @aliases NoneMultipliers
#' @exportClass NoneMultipliers
#'
#' @encoding latin1
#'
#' @keywords S4-classes
#'
#' @seealso \code{\link{getMultipliers-NoneMultipliers}}
################################################################################

setClass(
    Class = "NoneMultipliers",
    contains = "BootMultipliers"
)

setMethod(
    f = "initialize",
    signature = "NoneMultipliers",
    definition = function(.Object, l, N) {
      
      .Object@l <- 1
      .Object@N <- N
      
      # Return object
      return(.Object)
    }
)

################################################################################
#' Get Multipliers for the None Multiplier.
#'
#' @name getMultipliers-NoneMultipliers
#' @aliases getMultipliers,NoneMultipliers-method
#' 
#' @param object a \code{NoneMultipliers} object; used to specify the
#'                parameters \code{N}, \code{l} and the type of the
#'                scheme to generate the dependent multipliers.
#' @param B Number of independent repetitions to bootstrap.
#'
#' @return a matrix of dimension \code{[N,B]} with all 1.
################################################################################

setMethod(f = "getMultipliers",
    signature = "NoneMultipliers",
    definition = function(object, B=1) {
      
      N <- object@N
      l <- object@l
      
      return( matrix(1, nrow=N, ncol=B) )
      
    }
)


################################################################################
#' Create an instance of the \code{NoneMultipliers} class.
#'
#' l will be chosen as 1 by default.
#' 
#' @name NoneMultipliers-constructor
#' @aliases noneMultipliers
#' @export
#'
#' @keywords Constructors
#'
#' @param N   length of the sequence of dependent multipliers
#'
#' @return Returns an instance of \code{NoneMultipliers}.
################################################################################

noneMultipliers <- function( N ) {
  
  # TODO: Add more checks
  if (!(is.wholenumber(N) && N > 0)) {
    stop("'N' needs to be specified as integers with N > 0")
  }
  
  obj <- new(
      Class = "NoneMultipliers",
      l = 1,
      N = N
  )
  
  return(obj)
}
