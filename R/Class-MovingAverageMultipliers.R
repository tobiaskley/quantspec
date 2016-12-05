#' @include Class-BootMultipliers.R
NULL

################################################################################
#' Class for Moving Average Multiplier implementation.
#'
#' \code{MovingAverageMultipliers} is an S4 class that implements the moving blocks
#' bootstrap described in B{\"u}cher and Kojadinovic (2016), Section 5.2.1.
#'
#' \code{MovingAverageMultipliers} extends the S4 class
#' \code{\link{BootMultipliers}} and the remarks made in its documentation
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
        kappa = "function",
        distrInnov = "function"
    ),
    contains = "BootMultipliers"
)

setMethod(
    f = "initialize",
    signature = "MovingAverageMultipliers",
    definition = function(.Object, l, N, kappa, distrInnov) {
      
      .Object@l <- l
      .Object@N <- N
      .Object@kappa <- kappa
      .Object@distrInnov <- distrInnov
      
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
      distrInnov <- object@distrInnov
      
      b <- (l+1)/2
      w <- kappa( (1:l - b)/b )
      
      Z <- matrix( distrInnov( B * (N + 2*b - 2) ), ncol=B)
      
      #Z2 <- Z[,2]
      
      res  <- c()
      for (i in 1:N) {
        # res <- c(res, sum(w * Z2[i + 1:l - 1]))
        res <- rbind(res, w %*% Z[i + 1:l - 1, ])
      }
      
      return(res / sum(w^2)^(1/2) )
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

movingAverageMultipliers <- function( l, N, kappa, distrInnov = rnorm ) {
  
  # TODO: Add more checks
  if (!(is.wholenumber(l) && is.wholenumber(N) && 0 < l && l <= N)) {
    stop("'l' and 'N' need to be specified as integers with 0 < l <= N")
  }
  
  # enlarge l by +1, if even
  if (l %% 2 == 0) {
    message("'l' enlarged by +1, as l needs to be odd.")
    l <- l + 1
  }
  
  obj <- new(
      Class = "MovingAverageMultipliers",
      l = l,
      N = N,
      kappa = kappa,
      distrInnov = distrInnov
  )
  
  return(obj)
}
