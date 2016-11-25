#' @include Class-DependentMultipliers.R
NULL

################################################################################
#' Class for Covariance Matrix Multiplier implementation.
#'
#' \code{CovarianceMatrixMultipliers} is an S4 class that implements the moving blocks
#' bootstrap described in B{\"u}cher and Kojadinovic (2016), Section 5.2.2.
#'
#' \code{CovarianceMatrixMultipliers} extends the S4 class
#' \code{\link{DependentMultipliers}} and the remarks made in its documentation
#' apply here as well.
#'
#' TODO: Add more description here!
#'
#' @name   CovarianceMatrixMultipliers-class
#' @aliases CovarianceMatrixMultipliers
#' @exportClass CovarianceMatrixMultipliers
#'
#' @encoding latin1
#'
#' @keywords S4-classes
#' 
#' @slot phi function from assumption (M3)
#'
#' @seealso \code{\link{getMultipliers-CovarianceMatrixMultipliers}}
#'
#' @references
#' B{\"u}cher, A., & Kojadinovic, I. (2016).
#' A dependent multiplier bootstrap for the sequential empirical copula
#' process under strong mixing. \emph{Bernoulli}, \bold{22}(2), 927--968.
################################################################################

setClass(
    Class = "CovarianceMatrixMultipliers",
    representation=representation(
        phi = "function"
    ),
    contains = "DependentMultipliers"
)

setMethod(
    f = "initialize",
    signature = "CovarianceMatrixMultipliers",
    definition = function(.Object, l, N, phi) {
      
      .Object@l <- l
      .Object@N <- N
      .Object@N <- phi
      
      # Return object
      return(.Object)
    }
)

################################################################################
#' Get Multipliers for the Covariance Matrix Multiplier.
#'
#' @name getMultipliers-CovarianceMatrixMultipliers
#' @aliases getMultipliers,CovarianceMatrixMultipliers-method
#' 
#' @importFrom stats runif
#'
#' @param object a \code{CovarianceMatrixMultipliers} object; used to specify the
#'                parameters \code{N}, \code{l}, \code{phi} and the type of the
#'                scheme to generate the dependent multipliers.
#' @param B Number of independent repetitions to bootstrap.
#'
#' @return a matrix of dimension \code{[N,B]} where each column gives a sequence
#'         of multipliers generated according to the moving average scheme.
################################################################################

setMethod(f = "getMultipliers",
    signature = "CovarianceMatrixMultipliers",
    definition = function(object, B=1) {
      
      N <- object@N
      l <- object@l
      phi <- object@phi
#      nBlocks <- ceiling(N/l)
#      
#      positions <- c()
#      
#      for (b in 1:B) {
#        blocks <- matrix(ncol=nBlocks, nrow=l)
#        blocks[1,] <- floor(runif(n=nBlocks, min=1,max=N-l+1))
#        if (l > 1) {
#          for (i in 2:l) {
#            blocks[i,] <- blocks[1,]+i-1
#          }
#        }
#        positions <- c(positions,as.vector(blocks)[1:N])
#      }
#      
#      return(matrix(positions,nrow=N))
    }
)


################################################################################
#' Create an instance of the \code{\link{CovarianceMatrixMultipliers}} class.
#'
#' @name CovarianceMatrixMultipliers-constructor
#' @aliases covarianceMatrixMultipliers
#' @export
#'
#' @keywords Constructors
#'
#' @param l   the range of dependence of the multipliers (they are l-dependent)
#' @param N   length of the sequence of dependent multipliers
#' @param phi function from assumption (M3)
#'
#' @return Returns an instance of \code{CovarianceMatrixMultipliers}.
################################################################################

covarianceMatrixMultipliers <- function( l, N, phi ) {
  
  # TODO: Add more checks
  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }
  
  obj <- new(
      Class = "CovarianceMatrixMultipliers",
      l = l,
      N = N,
      phi = phi
  )
  
  return(obj)
}
