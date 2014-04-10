\docType{class}
\name{SpecDistrWeight-class}
\alias{SpecDistrWeight}
\alias{SpecDistrWeight-class}
\title{Class for weights to estimate integrated spectral density kernels.}
\description{
\code{SpecDistrWeight} is an S4 class that implements a
weighting function given by \deqn{W_n(\alpha) := I\{\alpha
\leq 0\}}{Wn(a) := I{a <= 0}}.
}
\details{
At position \code{k} the value \eqn{W_n(2\pi
(k-1)/n}{Wn(2pi(k-1)/n)} is stored [in a vector
\code{values} nested inside \code{env}] for
\code{k=1,\ldots,T}. The number \code{length(values)} of
Fourier frequencies for which \eqn{W_n}{Wn} will be
evaluated may be set on construction or updated when
evoking the method \code{getValues}.
}
\keyword{S4-classes}

