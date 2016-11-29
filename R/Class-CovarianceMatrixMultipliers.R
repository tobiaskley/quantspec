#' @include Class-BootMultipliers.R
NULL

################################################################################
#' Class for Covariance Matrix Multiplier implementation.
#'
#' \code{CovarianceMatrixMultipliers} is an S4 class that implements the moving blocks
#' bootstrap described in B{\"u}cher and Kojadinovic (2016), Section 5.2.2.
#'
#' \code{CovarianceMatrixMultipliers} extends the S4 class
#' \code{\link{BootMultipliers}} and the remarks made in its documentation
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
#' @slot distrInnov function with one argument n that generates the independent
#'                  innovations
#' @slot sqrt_Sigma the positive definite square root of Sigma_n
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
        phi = "function",
        distrInnov = "function",
        sqrt_Sigma = "matrix"
    ),
    contains = "BootMultipliers"
)

#' @importFrom stats toeplitz
setMethod(
    f = "initialize",
    signature = "CovarianceMatrixMultipliers",
    definition = function(.Object, l, N, phi, distrInnov) {
      
      .Object@l <- l
      .Object@N <- N
      .Object@phi <- phi
      .Object@distrInnov <- distrInnov
      
      Sigma <- toeplitz(phi( (0:(N-1))/l ))
      
      ## find sqrtm of Sigma
      ## Version 1
      e <- eigen(Sigma)
      V <- e$vectors
      
      .Object@sqrt_Sigma <- V %*% diag(sqrt(e$values)) %*% t(V)
      
#      ## Version 2
#      clChol <- chol(Sigma)
#      svdCC <- svd(t(clChol))
#      
#      sqrt_Sigma2 <- svdCC$u %*% diag(svdCC$d) %*% t(svdCC$u)

#      ## Version 3
#      require(expm)
#      sqrt_Sigma3 <- sqrtm(Sigma)
      
      
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
      #l <- object@l
      #phi <- object@phi
      distrInnov <- object@distrInnov
      
      sqrt_Sigma <- object@sqrt_Sigma
      
      return( sqrt_Sigma %*% matrix(distrInnov(N*B), ncol=B) )
      
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
#' @param distrInnov function with one argument n that generates the independent
#'                   innovations
#'
#' @return Returns an instance of \code{CovarianceMatrixMultipliers}.
################################################################################

covarianceMatrixMultipliers <- function( l, N, phi, distrInnov ) {
  
  # TODO: Add more checks
  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }
  
  obj <- new(
      Class = "CovarianceMatrixMultipliers",
      l = l,
      N = N,
      phi = phi,
      distrInnov = distrInnov
  )
  
  return(obj)
}
