\name{SmoothedPG-constructors}
\alias{SmoothedPG-constructors}
\alias{smoothedPG}
\title{Create an instance of the \code{SmoothedPG} class.}
\usage{
smoothedPG(object, frequencies = 2 * pi/length(object) * 0:(length(object) -
  1), levels.1 = 0.5, levels.2 = levels.1, isRankBased = TRUE,
  type = c("clipped", "qr"), type.boot = c("none", "mbb"),
  parallel = FALSE, B = 0, l = 1, weight = kernelWeight())
}
\arguments{
  \item{object}{a time series (\code{numeric}, \code{ts},
  or \code{zoo} object) from which to determine the
  smoothed periodogram; alternatively a
  \code{\link{QuantilePG}} object can be supplied.}

  \item{isRankBased}{If true the time series is first
  transformed to pseudo data.}

  \item{levels.1}{A vector of length \code{K1} containing
  the levels \code{x1} at which the SmoothedPG is to be
  determined.}

  \item{levels.2}{A vector of length \code{K2} containing
  the levels \code{x2}.}

  \item{frequencies}{A vector containing frequencies at
  which to determine the smoothed periodogram.}

  \item{type}{A flag to choose the type of the estimator.
  Can be either \code{"clipped"} or \code{"qr"}. In the
  first case \code{\link{ClippedFT}} is used as a frequency
  representation, in the second case
  \code{\link{QRegEstimator}} is used.}

  \item{B}{number of bootstrap replications}

  \item{l}{(expected) length of blocks}

  \item{type.boot}{A flag to choose a method for the block
  bootstrap; currently two options are implemented:
  \code{"none"} and \code{"mbb"} which means to do a moving
  blocks bootstrap with \code{B} and \code{l} as
  specified.}

  \item{parallel}{a flag to allow performing parallel
  computations, where possible.}

  \item{weight}{Object of type \code{\link{Weight}} to be
  used for smoothing.}
}
\value{
Returns an instance of \code{SmoothedPG}.
}
\description{
A \code{SmoothedPG} object can be created from either
\itemize{ \item a \code{numeric}, a \code{ts}, or a
\code{zoo} object \item a \code{QuantilePG} object. } If a
\code{QuantilePG} object is used for smoothing, only the
\code{weight}, \code{frequencies} and \code{levels.1} and
\code{levels.2} parameters are used; all others are
ignored. In this case the default values for the levels are
the levels of the \code{QuantilePG} used for smoothing. Any
subset of the levels available there can be chosen.
}
\details{
The parameter \code{type.boot} can be set to choose a block
bootstrapping procedure. If \code{"none"} is chosen, a
moving blocks bootstrap with \code{l=length(Y)} and
\code{N=length(Y)} would be done. Note that in that case
one would also chose \code{B=0} which means that
\code{getPositions} would never be called. If \code{B>0}
then each bootstrap replication would be the undisturbed
time series.
}
\examples{
Y <- rnorm(16)
levels.1 <- c(0.25,0.5,0.75)
weight <- kernelWeight(W=W0)

# Version 1a of the constructor -- for numerics:
sPG.ft <- smoothedPG(Y, levels.1 = levels.1, weight = weight, type="clipped")
sPG.qr <- smoothedPG(Y, levels.1 = levels.1, weight = weight, type="qr")

# Version 1b of the constructor -- for ts objects:
sPG.ft <- smoothedPG(wheatprices, levels.1 = c(0.05,0.5,0.95), weight = weight)

# Version 1c of the constructor -- for zoo objects:
sPG.ft <- smoothedPG(sp500, levels.1 = c(0.05,0.5,0.95), weight = weight)

# Version 2 of the constructor:
qPG.ft <- quantilePG(Y, levels.1 = levels.1, type="clipped")
sPG.ft <- smoothedPG(qPG.ft, weight = weight)
qPG.qr <- quantilePG(Y, levels.1 = levels.1, type="qr")
sPG.qr <- smoothedPG(qPG.qr, weight = weight)
}
\keyword{Constructors}

