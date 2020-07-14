context("ClippedFT")

test_that("clippedFT works as expected",{
      
      set.seed(2581)
  
      n <- 65   
  
      Y1 <- rnorm(n)
      Y2 <- rank(Y1)/n
      
      freq <- 2*pi*(0:(n-1))/n
      levels <- c(0.25,0.5)
      
      cFT1 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=TRUE)
      cFT2 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=FALSE)
      
      # Compute cFT 'by hand':
      res1 <- matrix(0, nrow = n, ncol = 2)
      res2 <- matrix(0, nrow = n, ncol = 2)
      for (l in 1:length(levels)) {
        for (f in 1:length(freq)) {
          ii <- complex(real = 0, imaginary = 1)
          res1[f,l] <- sum((Y2 <= levels[l]) * exp(-1*ii*(0:(n-1))*freq[f]))
          res2[f,l] <- sum((Y1 <= levels[l]) * exp(-1*ii*(0:(n-1))*freq[f]))
        }
      }
      
      V1 <- getValues(cFT1, frequencies=freq, levels=levels)
      V2 <- getValues(cFT2, frequencies=freq, levels=levels)
      
      expect_that(dim(V1),equals(c(n, 2, 1)))
      expect_that(dim(V2),equals(c(n, 2, 1)))
      expect_that(V1[,,1],equals(res1))
      expect_that(V2[,,1],equals(res2))
      
      #sum((V1[,,1] - res1)^2)
    }
)
