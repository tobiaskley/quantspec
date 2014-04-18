\name{ts-models}
\alias{ts-models}
\alias{ts1}
\alias{ts2}
\alias{ts3}
\title{Functions to simulate from the time series models in Kley et. al (2014).}
\usage{
ts1(n)

ts2(n)

ts3(n)
}
\arguments{
  \item{n}{length of the time series to be returned}
}
\description{
Functions to simulate from the time series models in Kley
et. al (2014).
}
\details{
\code{ts1} QAR(1) model from Dette et. al (2014+).

\code{ts2} AR(2) model from Li (2012):

\code{ts3} ARCH(1) model from Lee and Subba Rao (2012):
}
\examples{
# Plot sample paths:
plot(ts1(100), type="l")
plot(ts2(100), type="l")
plot(ts3(100), type="l")
}
\references{
Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2014+). Of
Copulas, Quantiles, Ranks and Spectra: an
\eqn{L_1}{L1}-approach to spectral analysis.
\emph{Bernoulli}, \bold{forthcoming}.

Li, T.-H. (2012). Quantile Periodograms. \emph{Journal of
the American Statistical Association}, \bold{107},
765--776.

Lee, J., & Subba Rao, S. (2012). The Quantile Spectral
Density and Comparison based Tests for Nonlinear Time
Series. \url{http://arxiv.org/abs/1112.2759}.
}

