\docType{methods}
\name{getStdError-QuantileSD}
\alias{getStdError,QuantileSD-method}
\alias{getStdError-QuantileSD}
\title{Get \code{stdError} from a quantile spectral density kernel}
\usage{
\S4method{getStdError}{QuantileSD}(object, frequencies = 2 * pi * (0:(object@N
  - 1))/object@N, levels.1 = getLevels(object, 1),
  levels.2 = getLevels(object, 2))
}
\arguments{
  \item{object}{\code{QuantileSD} of which to get the
  \code{stdError}}

  \item{frequencies}{a vector of frequencies for which to
  get the \code{stdError}}

  \item{levels.1}{the first vector of levels for which to
  get the \code{stdError}}

  \item{levels.2}{the second vector of levels for which to
  get the \code{stdError}}
}
\value{
Returns the array \code{stdError} that's a slot of
\code{object}.
}
\description{
The selection mechanism for frequencies and levels operates
in the same way as described in
\code{\link{getValues-QuantileSD}}. The format of the
output is also described there.
}
\keyword{Access-functions}

