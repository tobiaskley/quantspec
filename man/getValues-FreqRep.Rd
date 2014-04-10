\docType{methods}
\name{getValues-FreqRep}
\alias{getValues,FreqRep-method}
\alias{getValues-FreqRep}
\title{Get values from a frequency representation.}
\usage{
\S4method{getValues}{FreqRep}(object, frequencies = 2 * pi *
  (0:(length(object@Y) - 1))/length(object@Y), levels = object@levels)
}
\arguments{
  \item{object}{\code{FreqRep} of which to get the values}

  \item{frequencies}{a vector of frequencies for which to
  get the values}

  \item{levels}{a vector of levels for which to get the
  values}
}
\value{
Returns data from the array \code{values} that's a slot of
\code{object}.
}
\description{
For two vectors \code{frequencies} and \code{levels} the
values from an \code{object} of type \code{FreqRep} are
returned.
}
\details{
The two parameters \code{frequencies} and \code{levels} are
expected to be vectors of reals; an error is thrown
otherwise. If any of the \code{frequencies} or
\code{levels} requested is not available from \code{object}
a warning is issued, and the values with frequencies and
levels closest to the ones requested are returned. Note
that the frequencies are transformed to
\eqn{[0,\pi]}{[0,pi]} using
\code{\link{frequenciesValidator}} when checking if they
are available in \code{object}.

The returned array of \code{values} is of dimension
\code{[J,K,B+1]}, where \code{J=length(frequencies)},
\code{K=length(levels)}, and \code{B} denotes the value
stored in slot \code{B} of \code{object}. At position
\code{(j,k,b)} the returned value is the one corresponding
to \code{frequencies[j]} and \code{levels[k]} that are
closest to the \code{frequencies} and \code{levels}
available in \code{object}; \code{\link{closest.pos}} is
used to determine what closest to means.
}
\examples{
Y        <- rnorm(32)
freq     <- 2*pi*c(0:31)/32
levels   <- c(0.25,0.5,0.75)
cFT      <- clippedFT(Y, freq, levels)
V.all    <- getValues(cFT)
V.coarse <- getValues(cFT, frequencies = 2*pi*c(0:15)/16, levels = levels)
V.fine   <- getValues(cFT, frequencies = 2*pi*c(0:63)/64, levels = levels)
V.part   <- getValues(cFT, frequencies = 2*pi*c(0:16)/32, levels = c(0.25))
}
\keyword{Access-functions}

