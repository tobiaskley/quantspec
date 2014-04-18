\docType{class}
\name{QSpecQuantity-class}
\alias{QSpecQuantity}
\alias{QSpecQuantity-class}
\title{Class for a Quantile Spectral Estimator.}
\description{
\code{QSpecQuantity} is an S4 class that provides a common
interface to objects that are of the functional form
\eqn{f(\omega; x_1, x_2)}{f(w; x1, x2)}, where
\eqn{\omega}{w} is a frequency parameter and \eqn{x_1,
x_2}{x1, x2} are level parameters. For each combination of
parameters a complex number can be stored. Examples for
objects of this kind currently include the quantile (i. e.,
Laplace or copula) spectral density kernel [cf.
\code{\link{QuantileSD}} for an implementation], an
integrated version of the quantile spectral density kernels
[cf. \code{\link{IntegrQuantileSD}} for an implementation],
and estimators of it [cf. \code{\link{QuantilePG}} and
\code{\link{SmoothedPG}} for implementations].
}
\section{Slots}{

\describe{
\item{\code{values}}{The array holding the values \eqn{f(\omega; x_1, x_2)}{f(w; x1, x2)}.}

\item{\code{frequencies}}{The frequencies \eqn{\omega}{w} for which the values are
                   available.}

\item{\code{levels}}{A list of vectors containing the levels \eqn{x_i}{xi} serving
              as argument for the estimator.}
}}
\keyword{S4-classes}

