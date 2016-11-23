################################################################################
#' Simulation of an QAR(1) time series.
#'
#' Returns a simulated time series \eqn{(Y_t)} that fulfills the following equation:
#' \deqn{Y_t = \theta_1(U_t) Y_{t-1} + \theta_0(U_t),}
#' where \eqn{\theta_1} and \eqn{\theta_0} are parameters and \eqn{U_t} is
#' independent white noise with uniform \eqn{[0,1]} marginal distributions.
#'
#' @name ts-models-QAR1
#' @aliases QAR1
#' @export
#' 
#' @importFrom stats qnorm
#' @importFrom stats runif
#'
#' @param n length of the time series to be returned
#' @param th1 parameter function with one argument \code{u} defined on
#'            \eqn{[0,1]}
#' @param th0 parameter function with one argument \code{u} defined on
#'                \eqn{[0,1]}
#' @param overhead an integer specifying the ``warmup'' period to reach an
#'                  approximate stationary start for the times series
#' @return Returns an QAR(1) time series with specified parameters.
#'
#' @examples
#' plot(QAR1(100), type="l")
################################################################################
QAR1 <- function(n,th1=function(u){1.9*((u-0.5))}, overhead=1000, th0=qnorm) {
  Y <- rep(0,n+overhead)
  Y[1] <- th0(runif(1))
  for (t in 2:(n+overhead)) {
    U <- runif(1)
    Y[t] <- th1(U)*Y[t-1]+th0(U)
  }
  Y[(overhead+1):(overhead+n)]
}

################################################################################
#' Simulation of an AR(1) time series.
#'
#' Returns a simulated time series \eqn{(Y_t)} that fulfills the following equation:
#' \deqn{Y_t = a Y_{t-1} + \epsilon_t,}
#' where \eqn{a} is a parameter and \eqn{\epsilon_t} is independent white
#' noise with marginal distribution specified by the parameter \code{innov}.
#'
#' @name ts-models-AR1
#' @aliases AR1
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @param n length of the time series to be returned
#' @param a parameter of the model
#' @param overhead an integer specifying the ``warmup'' period to reach an
#'                  approximate stationary start for the times series
#' @param innov a function that generates a random number each time
#'               \code{innov(1)} is called; used to specify the distribution of
#'               the innovations; \code{rnorm} by default
#'
#' @return Returns an AR(1) time series with specified parameters.
#'
#' @examples
#' plot(AR1(100, a=-0.7), type="l")
################################################################################
AR1 <- function(n,a,overhead=500,innov=rnorm) {
  Y <- rep(0,n+overhead)
  Y[1] <- innov(1)
  for (t in 2:(n+overhead)) {
    Y[t] <- a*Y[t-1]+innov(1)
  }
  Y[(overhead+1):(overhead+n)]
}

################################################################################
#' Simulation of an AR(2) time series.
#'
#' Returns a simulated time series \eqn{(Y_t)} that fulfills the following equation:
#' \deqn{Y_t = a_1 Y_{t-1} + a_2 Y_{t-2} + \epsilon_t,}
#' where \eqn{a_1} and \eqn{a_2} are parameters and \eqn{\epsilon_t} is
#' independent white noise with marginal distribution specified by the
#' parameter \code{innov}.
#'
#' @name ts-models-AR2
#' @aliases AR2
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @param n length of the time series to be returned
#' @param a1 parameter
#' @param a2 parameter
#' @param overhead an integer specifying the ``warmup'' period to reach an
#'                  approximate stationary start for the times series
#' @param innov a function with one parameter \code{n} that yields \code{n}
#'               independent pseudo random numbers each time it is called.
#' @return Return an AR(2) time series with specified parameters.
#'
#' @examples
#' plot(AR2(100, a1=0, a2=0.5), type="l")
################################################################################
AR2 <- function(n,a1,a2,overhead=500,innov=rnorm) {
  Y <- rep(0,n+overhead)
  Y[1] <- innov(1)
  Y[2] <- innov(1)
  for (t in 3:(n+overhead)) {
    Y[t] <- a1*Y[t-1]+a2*Y[t-2]+innov(1)
  }
  Y[(overhead+1):(overhead+n)]
}

################################################################################
#' Simulation of an ARCH(1) time series.
#'
#' Returns a simulated time series \eqn{(Y_t)} that fulfills the following equation:
#' \deqn{Y_t = Z_t \sigma_t, \quad \sigma_t^2 = a_0 + a_1 Y_{t-1}^2 + \epsilon_t}
#' where \eqn{a_0} and \eqn{a_1} are parameters and \eqn{\epsilon_t} is
#' independent white noise with marginal distribution specified by the
#' parameter \code{innov}.
#'
#' @name ts-models-ARCH1
#' @aliases ARCH1
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @param n length of the time series to be returned
#' @param a0 parameter
#' @param a1 parameter
#' @param overhead an integer specifying the ``warmup'' period to reach an
#'                  approximate stationary start for the times series
#' @param innov a function with one parameter \code{n} that yields \code{n}
#'               independent pseudo random numbers each time it is called.
#' @return Return an ARCH(1) time series with specified parameters.
#'
#' @examples
#' plot(ARCH1(100, a0=1/1.9, a1=0.9), type="l")
#'
################################################################################
ARCH1 <- function(n,a0,a1,overhead=500,innov = rnorm) {
  Y <- rep(0,n+overhead)
  Y[1] <- 1
  for (t in 2:(n+overhead)) {
    Y[t] <- sqrt(a0 + a1*Y[t-1]^2) * innov(1)
  }
  Y[(overhead+1):(overhead+n)]
}

#ts1 <- function(n){AR1(n,-0.5,innov=function(n){rnorm(n)})}
#ts2 <- function(n){AR1(n,-0.3,innov=function(n){rt(n,1)})}


################################################################################
#' Functions to simulate from the time series models in Kley et. al (2016).
#'
#' @details
#' \code{ts1} QAR(1) model from Dette et. al (2015).
#'
#' @name ts-models
#' @aliases ts1
#' @export
#' 
#' @importFrom stats qnorm
#'
#' @param n length of the time series to be returned
#'
#' @references
#' Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2015).
#' Of Copulas, Quantiles, Ranks and Spectra: an \eqn{L_1}{L1}-approach to
#' spectral analysis. \emph{Bernoulli}, \bold{21}(2), 781--831.
#' [cf. \url{http://arxiv.org/abs/1111.7205}]
#'
#' @examples
#' # Plot sample paths:
#' plot(ts1(100), type="l")
################################################################################
ts1 <- function(n){QAR1(n,th0=function(x){0.1*qnorm(x)})}

################################################################################
#' @details
#' \code{ts2} AR(2) model from Li (2012):
#'
#' @name ts-models
#' @aliases ts2
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @references
#' Li, T.-H. (2012). Quantile Periodograms.
#' \emph{Journal of the American Statistical Association}, \bold{107}, 765--776.
#'
#' @examples
#' plot(ts2(100), type="l")
################################################################################
ts2 <- function(n){
  #  r <- 0.6
  #  wc <- 2*pi*0.25
  #  a1 <- 2*r*cos(wc)
  #  a2 <- -1*r^2
  a1 <- 0
  a2 <- -0.36
  AR2(n,a1,a2,innov=function(n){rnorm(n)})
}


################################################################################
#' @details
#' \code{ts3} ARCH(1) model from Lee and Subba Rao (2012):
#'
#' @name ts-models
#' @aliases ts3
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @references
#' Lee, J., & Subba Rao, S. (2012).
#' The Quantile Spectral Density and Comparison based Tests for Nonlinear Time
#' Series. \url{http://arxiv.org/abs/1112.2759}.
#'
#' @examples
#' plot(ts3(100), type="l")
################################################################################

# Example from Subba Rao's paper
ts3 <- function(n){ARCH1(n,1/1.9,0.9,innov=function(n){rnorm(n)})}


################################################################################
#' Functions to simulate multivariate time series models from B{\"u}cher and
#' Kojadinovic (2016)
#' 
#' @details
#' In these bivariate time series models the cross sectional dependence is
#' introduced via a cross sectional dependence of the innovations. More precisely,
#' \eqn{(\varepsilon_{i1}, \varepsilon_{i2}) := (\Phi^{-1}(U_{i1}), \Phi^{-1}(U_{i2}))},
#' with \eqn{\Phi} the c.d.f. of the standard normal distribution and
#' \eqn{(U_{i1}, U_{i2})} bivariate i.i.d. from either a Clayton copula with
#' parameters \eqn{1} or \eqn{4} or Gumbel copula \eqn{1.5} or \eqn{3},
#' respectively.
#' 
#' The four models then have the following mechanisms to introduce serial
#' dependence: 
#' 
#' \code{mvts1} bivariate AR(1) model:
#' \deqn{X_{ij} = 0.5 X_{i-1,j} + \varepsilon_{ij},}
#'
#' @name ts-models-mv
#' @aliases mvts1
#' @export
#' 
#' @importFrom copula mvdc
#' @importFrom copula archmCopula
#' @importFrom copula rMvdc
#' @importFrom zoo zoo
#' 
#' @param n     length of the time series to be returned
#' @param innov flag to indicate whether the independent innovations are
#'              distributed according to Clayton copula with parameters 1 or 4
#'              or according to Gumbel copula with parameters 1.5 or 3.
#'              
#' @references
#' B{\"u}cher, A., & Kojadinovic, I. (2016).
#' A dependent multiplier bootstrap for the sequential empirical copula
#' process under strong mixing. \emph{Bernoulli}, \bold{22}(2), 927--968.
#'
#' @examples
#' require("zoo")
#' Y <- zoo(mvts1(100, innov = "Clayton4"))
#' plot(Y)
################################################################################

# Example from Axel's paper
mvts1 <- function(n = 100, innov = c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3")) {
  overhead <- 100
  innov <- match.arg(innov, c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3"))[1]
  switch(innov,
      "Clayton1" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 1),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Clayton4" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 4),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel15" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 1.5),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel3" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 3),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))}
  )
  eps <- rMvdc(n+overhead+1, innovDist)
  Y <- matrix(0, nrow = n+overhead, ncol = 2)
  Y[1, ] <- eps[1, ]
  for (t in 2:(n+overhead)) {
    Y[t, ] <- 0.5 * Y[t-1, ] + eps[t, ]
  }
  Y[(overhead+1):(overhead+n), ]
}


################################################################################
#' @details
#' \code{mvts2} bivariate NAR(1) model:
#' \deqn{X_{ij} = 0.6 sin(X_{i-1,j}) + \varepsilon_{ij},} 
#'
#' @name ts-models-mv
#' @aliases mvts2
#' @export
#' 
#' @importFrom copula mvdc
#' @importFrom copula archmCopula
#' @importFrom copula rMvdc
#' @importFrom zoo zoo
#' 
#' @examples
#' require("zoo")
#' Y <- zoo(mvts2(100, innov = "Gumbel15"))
#' plot(Y)
################################################################################

# Example from Axel's paper
mvts2 <- function(n = 100, innov = c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3")) {
  overhead <- 100
  innov <- match.arg(innov, c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3"))[1]
  switch(innov,
      "Clayton1" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 1),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Clayton4" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 4),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel15" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 1.5),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel3" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 3),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))}
  )
  eps <- rMvdc(n+overhead+1, innovDist)
  Y <- matrix(0, nrow = n+overhead, ncol = 2)
  Y[1, ] <- eps[1, ]
  for (t in 2:(n+overhead)) {
    Y[t, ] <- 0.6 * sin(Y[t-1, ]) + eps[t, ]
  }
  Y[(overhead+1):(overhead+n), ]
}


################################################################################
#' @details
#' \code{mvts3} bivariate EXPAR(1) model:
#' \deqn{X_{ij} = (0.8 - 1.1 exp(-50 X^2_{i-1,j}) X_{i-1,j} + 0.1 \varepsilon_{ij},} 
#'
#' @name ts-models-mv
#' @aliases mvts3
#' @export
#' 
#' @importFrom copula mvdc
#' @importFrom copula archmCopula
#' @importFrom copula rMvdc
#' @importFrom zoo zoo
#' 
#' @examples
#' require("zoo")
#' Y <- zoo(mvts3(100, innov = "Gumbel15"))
#' plot(Y)
################################################################################

# Example from Axel's paper
mvts3 <- function(n = 100, innov = c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3")) {
  overhead <- 100
  innov <- match.arg(innov, c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3"))[1]
  switch(innov,
      "Clayton1" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 1),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Clayton4" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 4),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel15" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 1.5),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel3" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 3),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))}
  )
  eps <- rMvdc(n+overhead+1, innovDist)
  Y <- matrix(0, nrow = n+overhead, ncol = 2)
  Y[1, ] <- eps[1, ]
  for (t in 2:(n+overhead)) {
    Y[t, ] <- (0.8 - 1.1 * exp(-50 * Y[t-1, ]^2)) * Y[t-1, ] + 0.1 * eps[t, ]
  }
  Y[(overhead+1):(overhead+n), ]
}

################################################################################
#' @details
#' \code{mvts4} bivariate GARCH(1,1) model:
#' \deqn{X_{ij} = \sigma_{ij} \varepsilon_{ij}}
#' \deqn{\varepsilon_{i1}^2 = 0.012 + 0.919 \varepsilon_{i-1,1}^2 + 0.072 \varepsilon_{i-1,1}^2,}
#' #' and
#' \deqn{\varepsilon_{i2}^2 = 0.037 + 0.868 \varepsilon_{i-1,2}^2 + 0.115 \varepsilon_{i-1,2}^2.}
#' 
#' @name ts-models-mv
#' @aliases mvts4
#' @export
#' 
#' @importFrom copula mvdc
#' @importFrom copula archmCopula
#' @importFrom copula rMvdc
#' @importFrom zoo zoo
#' 
#' @examples
#' require("zoo")
#' Y <- zoo(mvts4(100, innov = "Gumbel15"))
#' plot(Y)
################################################################################

# Example from Axel's paper
mvts4 <- function(n = 100, innov = c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3")) {
  overhead <- 100
  innov <- match.arg(innov, c("Clayton1", "Clayton4", "Gumbel15", "Gumbel3"))[1]
  switch(innov,
      "Clayton1" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 1),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Clayton4" = {
        innovDist <- mvdc(copula = archmCopula(family = "clayton", dim = 2, param = 4),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel15" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 1.5),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))},
      "Gumbel3" = {
        innovDist <- mvdc(copula = archmCopula(family = "gumbel", dim = 2, param = 3),
            margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1)))}
  )
  eps <- rMvdc(n+overhead+1, innovDist)
  Y <- matrix(0, nrow = n+overhead, ncol = 2)
  sigma1 <- sqrt( 0.012 / (1 - 0.919 - 0.072) )
  sigma2 <- sqrt( 0.037 / (1 - 0.868 - 0.115) )

  Y[1, 1] <- sigma1 * eps[1, 1]
  Y[1, 2] <- sigma2 * eps[1, 2]
  for (t in 2:(n+overhead)) {
    sigma1 <- 0.012 + 0.919 * sigma1 + 0.072 * eps[t-1,1]^2
    sigma2 <- 0.037 + 0.868 * sigma2 + 0.115 * eps[t-1,2]^2
    Y[t, 1] <- sigma1 * eps[t, 1]
    Y[t, 2] <- sigma2 * eps[t, 2]
  }
  Y[(overhead+1):(overhead+n), ]
}