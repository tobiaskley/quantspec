#' @include generics.R
NULL

################################################################################
#' Class for Generation of Multipliers to be used for bootstrapping.
#'
#' \code{DependentMultipliers} is an S4 class that provides a common interface
#' to different algorithms that can be used for implementation of a dependent
#' multiplier sequence \eqn{\xi_{1,n}, \ldots, \xi_{n,n}}, satisfying the
#' conditions (M1)--(M3) stated in B{\"u}cher and Kojadinovic (2016), p. 930.
#'
#' After initialization the multipliers can be retrieved by applying
#' \code{getMultipliers} to the object.
#'
#' Different appraoches to generate the multipliers are implemented by creating
#' a subclass together with a \code{getMultipliers} method that contains the
#' implementation of dependent multiplier scheme.
#'
#' Currently the following implementations are available:
#'
#' \itemize{
#' 		\item \code{\link{MovingAverageMultipliers}}
#'           and \code{\link{getMultipliers-MovingAverageMultipliers}},
#' 		\item \code{\link{CovarianceMatrixMultipliers}}
#'           and \code{\link{getMultipliers-CovarianceMatrixMultipliers}},
#' }
#'
#' @name   BootMultipliers-class
#' @aliases BootMultipliers
#' @exportClass BootMultipliers
#'
#' @keywords S4-classes
#'
#' @slot l   the range of dependence of the multipliers (they are l-dependent)
#' @slot N   length of the sequence of dependent multipliers
#'
#' @references
#' B{\"u}cher, A., & Kojadinovic, I. (2016).
#' A dependent multiplier bootstrap for the sequential empirical copula
#' process under strong mixing. \emph{Bernoulli}, \bold{22}(2), 927--968.
################################################################################

setClass(
    Class = "BootMultipliers",
    representation=representation(
        l = "numeric",
        N = "numeric"
    )
)
