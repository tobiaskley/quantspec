################################################################################
# Simulation study, analyzing the QAR(1) model from Dette et. al (2014+):
# -----------------------------------------------------------------------
# In this demo a number of R quantile periodograms and smoothed quantile
# periodograms are computed from R independent, simulated time series from the
# QAR(1) model with parameters as in the above paper.
#
# The root integrated mean squared errors are stored and can be used for
# comparison of the estimators.
#
# Finally, plots for the last simulation run are shown.
################################################################################

ts <- ts3
type <- "copula"
levels <- c(0.25, 0.5, 0.75)
freq <- 2*pi*(1:32)/64
R <- 1

J <- length(freq)
K <- length(levels)

# First determine the copula spectral density
csd <- quantileSD(N=2^8, seed.init = 2581, type = type,
                  ts = ts, levels.1=levels, R = 100)

# init an array for the root integrated mean squared errors
# (1: CR periodogram, 2: rank-based Laplace periodogram,
#  3: smoothed CR periodogram, 4: smoothed rank-based Laplace pg.)
rimse  <- array(0, dim=c(4,J,K,K))

weight <- kernelWeight(W=W1, N=32, bw=0.3)

trueV <- getValues(csd, frequencies=freq)
for (i in 1:R) {
  Y <- ts3(64)
  
  CR <- quantilePG(Y, levels.1=levels, type="clipped")
  rimse[1,,,] <- rimse[1,,,] +
      abs(getValues(CR, frequencies=freq)[,,,1]-trueV)^2
      
  
  LP <- quantilePG(Y, levels.1=levels, type="qr")
  rimse[2,,,] <- rimse[2,,,] + 
      abs(getValues(LP, frequencies=freq)[,,,1]-trueV)^2
  
  sCR <- smoothedPG(CR, weight=weight) # , frequencies=freq
  rimse[3,,,] <- rimse[3,,,] + 
      abs(getValues(sCR, frequencies=freq)[,,,1]-trueV)^2
  
  sLP <- smoothedPG(LP, weight=weight) # , frequencies=freq
  rimse[4,,,] <- rimse[4,,,] + 
      abs(getValues(sLP, frequencies=freq)[,,,1]-trueV)^2
}

rimse <- sqrt(apply(rimse, c(1,3,4), mean) / R)

# Inspect the rimse, but note that to have reliable approximations of
# the true rimse the simulations should run much longer.
for (i in 1:4) {rimse[i,,]}

# Finally take a look a the last simulated periodograms and see that
# time series of length 64 are not quiet long enought to yield
# reliable estimators.

plot(sCR, qsd=csd, frequencies = freq)
plot(sCR, qsd=csd, plotPG=TRUE, frequencies = freq)
plot(sLP, qsd=csd, frequencies = freq,
    ptw.CIs = 0, type.scaling="real-imaginary")
plot(sLP, qsd=csd, plotPG=TRUE, frequencies = freq,
    ptw.CIs = 0, type.scaling="real-imaginary")
