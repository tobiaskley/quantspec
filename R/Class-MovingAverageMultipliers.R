#' @include Class-DependentMultipliers.R
NULL

################################################################################
#' Class for Moving Average Multiplier implementation.
#'
#' \code{MovingAverageMultipliers} is an S4 class that implements the moving blocks
#' bootstrap described in B{\"u}cher and Kojadinovic (2016), Section 5.2.1.
#'
#' \code{MovingAverageMultipliers} extends the S4 class
#' \code{\link{DependentMultipliers}} and the remarks made in its documentation
#' apply here as well.
#'
#' TODO: Add more description here!
#'
#' @name   MovingAverageMultipliers-class
#' @aliases MovingAverageMultipliers
#' @exportClass MovingAverageMultipliers
#'
#' @encoding latin1
#'
#' @keywords S4-classes
#' 
#' @slot kappa      function from Section 5.2.1 (B{\"u}cher and Kojadinovic, 2016)
#' @slot distrInnov function with one argument n that generates the independent
#'                  innovations
#'
#' @seealso \code{\link{getMultipliers-MovingAverageMultipliers}}
#'
#' @references
#' B{\"u}cher, A., & Kojadinovic, I. (2016).
#' A dependent multiplier bootstrap for the sequential empirical copula
#' process under strong mixing. \emph{Bernoulli}, \bold{22}(2), 927--968.
################################################################################

setClass(
    Class = "MovingAverageMultipliers",
    representation=representation(
        kappa = "function"
    ),
    contains = "DependentMultipliers"
)

setMethod(
    f = "initialize",
    signature = "MovingAverageMultipliers",
    definition = function(.Object, l, N, phi) {
      
      .Object@l <- l
      .Object@N <- N
      .Object@N <- phi
      
      # Return object
      return(.Object)
    }
)

################################################################################
#' Get Multipliers for the Moving Average Multiplier.
#'
#' @name getMultipliers-MovingAverageMultipliers
#' @aliases getMultipliers,MovingAverageMultipliers-method
#' 
#' @importFrom stats runif
#'
#' @param object a \code{MovingAverageMultipliers} object; used to specify the
#'                parameters \code{N}, \code{l}, \code{phi} and the type of the
#'                scheme to generate the dependent multipliers.
#' @param B Number of independent repetitions to bootstrap.
#'
#' @return a matrix of dimension \code{[N,B]} where each column gives a sequence
#'         of multipliers generated according to the moving average scheme.
################################################################################

setMethod(f = "getMultipliers",
    signature = "MovingAverageMultipliers",
    definition = function(object, B=1) {
      
      N <- object@N
      l <- object@l
      kappa <- object@kappa
      
      b <- (l+1)/2
      w <- kappa( (1:l - b)/b )
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
#' Create an instance of the \code{\link{MovingAverageMultipliers}} class.
#'
#' @name MovingAverageMultipliers-constructor
#' @aliases movingAverageMultipliers
#' @export
#'
#' @keywords Constructors
#'
#' @param l     the range of dependence of the multipliers (they are l-dependent);
#'              for this scheme l needs to be chosen as an odd number and >= 3
#' @param N     length of the sequence of dependent multipliers
#' @param kappa function from Section 5.2.1 (B{\"u}cher and Kojadinovic, 2016)
#' @param distrInnov function with one argument n that generates the independent
#'                   innovations
#'
#' @return Returns an instance of \code{MovingAverageMultipliers}.
################################################################################

movingAverageMultipliers <- function( l, N, kappa, innovDistr = rnorm ) {
  
  # TODO: Add more checks
  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }
  
  obj <- new(
      Class = "MovingAverageMultipliers",
      l = l,
      N = N,
      kappa = kappa,
      innovDistr = innovDistr
  )
  
  return(obj)
}
