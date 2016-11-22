context("ClippedFT")

test_that("clippedFT works as expected",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      #Y2 <- rank(Y1)/64
      Fn <- ecdf(Y1)
      Y2 <- Fn(Y1)
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25, 0.5)
      
      cFT1 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=TRUE)
      cFT2 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=FALSE)
      
      # Compute cFT 'by hand':
      res1 <- matrix(0, nrow=64, ncol=2)
      res2 <- matrix(0, nrow=64, ncol=2)
      for (l in 1:length(levels)) {
        for (f in 1:length(freq)) {
          ii <- complex(real = 0, imaginary = 1)
          res1[f,l] <- sum((Y2 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
          res2[f,l] <- sum((Y1 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
        }
      }
      
      V1 <- getValues(cFT1, frequencies=freq, levels=levels)
      V2 <- getValues(cFT2, frequencies=freq, levels=levels)
      
      expect_that(dim(V1),equals(c(64,2,1)))
      expect_that(dim(V2),equals(c(64,2,1)))
      expect_that(V1[,,1],equals(res1))
      expect_that(V2[,,1],equals(res2))
      
    }
)

test_that("clippedFT works as expected (in the presence of ties)",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      Y1[31:35] <- Y1[1:5] 
      #Y2 <- rank(Y1)/64
      Fn <- ecdf(Y1)
      Y2 <- Fn(Y1)
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25, 10.7/64)
      
      cFT1 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=TRUE)
      cFT2 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=FALSE)
      
      # Compute cFT 'by hand':
      res1 <- matrix(0, nrow=64, ncol=2)
      res2 <- matrix(0, nrow=64, ncol=2)
      for (l in 1:length(levels)) {
        for (f in 1:length(freq)) {
          ii <- complex(real = 0, imaginary = 1)
          res1[f,l] <- sum((Y2 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
          res2[f,l] <- sum((Y1 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
        }
      }
      
      V1 <- getValues(cFT1, frequencies=freq, levels=levels)
      V2 <- getValues(cFT2, frequencies=freq, levels=levels)
      
      expect_that(dim(V1),equals(c(64,2,1)))
      expect_that(dim(V2),equals(c(64,2,1)))
      expect_that(V1[,,1],equals(res1))
      expect_that(V2[,,1],equals(res2))
      
    }
)

test_that("clippedFT works as expected (with bootstrap)",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25, 10.7/64)
      
      set.seed(2581)
      cFT1 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=FALSE, type.boot = "mbb", B = 1, l = 8)
      
      #set.seed(2581)
      #cFT2 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=TRUE, type.boot = "mbb", B = 1, l = 8, resampleEcdf = FALSE)
      
      set.seed(2581)
      cFT3 <- clippedFT(Y1, frequencies = freq, levels = levels, isRankBased=TRUE, type.boot = "mbb", B = 1, l = 8) #, resampleEcdf = TRUE)
      
      set.seed(2581)
      pos <- getPositions(getBootPos(cFT3), 1)
      
      Fn <- ecdf(Y1)
      Y2 <- Fn(Y1[pos])
      
      Fn <- ecdf(Y1[pos])
      Y3 <- Fn(Y1[pos])
      
      # Compute cFT 'by hand':
      res1 <- matrix(0, nrow=64, ncol=2)
      res2 <- matrix(0, nrow=64, ncol=2)
      res3 <- matrix(0, nrow=64, ncol=2)
      for (l in 1:length(levels)) {
        for (f in 1:length(freq)) {
          ii <- complex(real = 0, imaginary = 1)
          res1[f,l] <- sum((Y1[pos] <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
          res2[f,l] <- sum((Y2 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
          res3[f,l] <- sum((Y3 <= levels[l]) * exp(-1*ii*(0:63)*freq[f]))
        }
      }
      
      V1 <- getValues(cFT1, frequencies=freq, levels=levels)
      #V2 <- getValues(cFT2, frequencies=freq, levels=levels)
      V3 <- getValues(cFT3, frequencies=freq, levels=levels)
      
      expect_that(dim(V1),equals(c(64, 2, 2)))
      #expect_that(dim(V2),equals(c(64, 2, 2)))
      expect_that(dim(V3),equals(c(64, 2, 2)))
      expect_that(V1[,,2],equals(res1))
      #expect_that(V2[,,2],equals(res2))
      expect_that(V3[,,2],equals(res3))
      
    }
)
