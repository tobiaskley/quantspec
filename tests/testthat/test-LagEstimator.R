context("LagEstimator")

test_that("lagEstimator works as expected",{
  
  set.seed(9247)
  Bn = 10
  Y <- rnorm(100)
  levels.1 <- c(0.1,0.5,0.9)
  weight <- lagKernelWeight(W = WParzen,  bw = Bn, K = length(Y))
  lagOp <- clippedCov(Y,levels.1 = levels.1)
  lagEst <- lagEstimator(lagOp,weight = weight)
  V <- getValues(lagEst)
  
  
  #Calculate lag-window estimator by hand:
  
  # weight function
  K_Parzen_aux <- function (u) {
    if (abs(u) <= 1) {
      if(abs(u)<=0.5){return(1-6*u^2+6*abs(u)^3)}    
      else{return(2*(1-abs(u))^3)}
      
    } else {
      return(0)
    }
  }
  
  Kernel <- Vectorize(K_Parzen_aux)
  
  
  K <- length(levels.1)
  n <- length(Y)
  Q <- quantile(Y, probs=levels.1)
  
  # Compute the matrix of indicators and the cross-correlations Fhat
  
  IndMatrix <- matrix(0, nrow=n, ncol=K)
  
  for (k in 1:K) {
    IndMatrix[,k] <- (Y <= Q[k]) - levels.1[k]
    #IndMatrix[,k] <- IndMatrix[,k] - mean(IndMatrix[,k]) 
  }
  
  Fhat <- acf(IndMatrix, type="covariance", lag.max = n-1, plot = FALSE, demean = FALSE)
  
  # Define and fill a result vector:
  
  res <- array(dim = c(n,K,K))
  
  for (k1 in 1:K) {
    for (k2 in k1:K) {
      A1 <- Fhat$acf[,k1,k2]
      A1 <- A1*Kernel((0:(n-1))/Bn)
      
      A2 <- Fhat$acf[,k2,k1]
      A2 <- A2*Kernel((0:(n-1))/Bn)
      
      A0 <- Fhat$acf[1,k2,k1]
      A0 <- A0*Kernel(0)
      
      res[,k1,k2] <- 1/(2*pi) * (fft(A1) + Conj(fft(A2)) - A0)
      if (k1 != k2) {res[,k2,k1] <- Conj(res[,k1,k2])}
    }
  }
  
  expect_equal(dim(V),c(100,3,3,1))
  expect_equal(V[,,,1],res)
}
)