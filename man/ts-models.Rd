\name{ts-models}
\alias{ts-models}
\alias{ts3}
\alias{ts5}
\alias{ts8}
\title{Functions to simulate from the time series models in Kley et. al (2014).}
\usage{
ts3(n)

ts5(n)

ts8(n)
}
\arguments{
  \item{n}{length of the time series to be returned}
}
\description{
Functions to simulate from the time series models in Kley
et. al (2014).
}
\details{
\code{ts3} QAR(1) model from Dette et. al (2014+).

\code{ts5} AR(2) model from Li (2012):

\code{ts8} ARCH(1) model from Lee and Subba Rao (2012):
}
\examples{
# Plot sample paths:
plot(ts3(100), type="l")
plot(ts5(100), type="l")
plot(ts8(100), type="l")
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

