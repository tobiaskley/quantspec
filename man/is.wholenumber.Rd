\name{is.wholenumber}
\alias{is.wholenumber}
\title{Checks whether \code{x} contains integer numbers.}
\usage{
is.wholenumber(x, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{x}{a vector to be checked for integers}

  \item{tol}{an optional parameter specifying to which
  precision the check is to be performed.}
}
\value{
Returns the time series as a numeric.
}
\description{
Borrowed from the example in \code{\link[base]{integer}}.
}
\examples{
\dontrun{
is.wholenumber(1) # is TRUE
(x <- seq(1, 5, by = 0.5) )
is.wholenumber( x ) #-->  TRUE FALSE TRUE ...
}
}
\keyword{internals}

