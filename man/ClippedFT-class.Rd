\docType{class}
\name{ClippedFT-class}
\alias{ClippedFT}
\alias{ClippedFT-class}
\title{Class for Fourier transform of the clipped time series.}
\description{
\code{ClippedFT} is an S4 class that implements the
necessary calculations to determine the Fourier transform
of the clipped time series. As a subclass to
\code{\link{FreqRep}} it inherits slots and methods defined
there; it servers as a frequency representation of a time
series as described in Kley et. al (2014+).
}
\details{
For each frequency \eqn{\omega} from \code{frequencies} and
level \code{q} from \code{levels} the statistic
\deqn{\sum_{t=0}^{n-1} I\{Y_t \leq q\} \mbox{e}^{-\mbox{i}
\omega t}} is determined and stored to the array
\code{values}. Internally the methods
\code{\link[stats]{mvfft}} and \code{\link[stats]{fft}} are
used to achieve good performance.

Note that, all remarks made in the documentation of the
super-class \code{\link{FreqRep}} apply.
}
\references{
Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2014).
Quantile Spectral Processes: Asymptotic Analysis and
Inference. \url{http://arxiv.org/abs/1401.8104}.
}
\seealso{
For an example see \code{\link{FreqRep}}.
}
\keyword{S4-classes}

