\docType{data}
\name{data-sp500}
\alias{data-sp500}
\alias{sp500}
\title{S&P 500: Standard and Poor's 500 stock index, 2004--2006}
\format{A univariate time series with 755 observations; a \code{zoo} object}
\description{
Contains the returns of the S&P 500 stock index for the
years 2004--2006. The returns were computed as
\code{(Adjusted.Close-Open)/Open}.
}
\details{
The data was downloaded from the Yahoo! Finance Website.
}
\examples{
plot(sp500)
}
\references{
Yahoo! Finance Website
\url{http://finance.yahoo.com/q/hp?s=^GSPC+Historical+Prices}
}
\keyword{data}

