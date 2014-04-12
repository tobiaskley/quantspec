## This script can be used to create and store a QuantileSD object

## Parameters for the simulation:
R <- 50                      # number of independent repetitions;
                             # R should be much larger than this in practice!
N <- 2^8                     # number of Fourier frequencies in [0,2pi)
ts <- ts1                    # time series model
levels <- seq(0.1,0.9,0.1)   # quantile levels
type <- "copula"             # copula, not Laplace, spectral density kernel
seed.init <- 2581            # seed for the pseudo random numbers

## Simulation takes place once the constructor is invoked
qsd <- quantileSD(N=N, seed.init = 2581, type = type,
    ts = ts, levels.1=levels, R = R)

## The simulated copula spectral density kernel can be called via
V1 <- getValues(qsd)

## It is also possible to fetch the result for only a few levels
levels.few <- c(0.2,0.5,0.7)
V2 <- getValues(qsd, levels.1=levels.few, levels.2=levels.few)

## If desired additional repetitions can be performed to yield a more precise
## simulation result by calling; here the number of independent runs is doubled.
sCSD <- increasePrecision(qsd,R)

## Often the result will be stored for later usage.  
save(qsd, file="QAR1.rdata")

## Take a brief look at the result of the simulation
plot(qsd, levels=levels.few)

## When plotting more than only few levels it may be a good idea to plot to
## another device; e. g., a pdf-file
K <- length(levels)
pdf("QAR1.pdf", width=2*K, height=2*K)
  plot(qsd)
dev.off()
