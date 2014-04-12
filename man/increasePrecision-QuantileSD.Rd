\docType{methods}
\name{increasePrecision-QuantileSD}
\alias{increasePrecision,QuantileSD-method}
\alias{increasePrecision-QuantileSD}
\title{Increase the precision of a \code{QuantileSD}}
\usage{
\S4method{increasePrecision}{QuantileSD}(object, R = 1, quiet = FALSE)
}
\arguments{
  \item{object}{The \code{\link{QuantileSD}} of which to
  increase the precision.}

  \item{R}{value of which to enlarge R}

  \item{quiet}{Don't report progress to console when
  computing the \code{R} independent quantile
  periodograms.}
}
\value{
Returns an \code{\link{QuantileSD}} object determined from
\code{oldR + R} independent repetitions.
}
\description{
The precision is increased by generating an additional
\code{R} \code{\link{QuantilePG}} objects (independent of
the previous ones) and then including them in the average.
}
\examples{
# First simulate a copula spectral density from R=20 independent runs.
csd <- quantileSD(N=2^9, ts=ts1, levels.1=c(0.25,0.5), type="copula", R=20)

# Check out the result:
getR(csd)
plot(csd)

# Now increase the number of independent simulation runs to 50.
csd <- increasePrecision(csd, R=30)

# Check out the (more precise) result:
getR(csd)
plot(csd)
}

