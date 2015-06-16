context("frequencyValidator")

test_that("frequencyValidator works as expected",{
			
			freq <- 2*pi*c(3,2,11,8,5)/10
			
			res <- frequenciesValidator(freq, N=10, steps=1:3)*10/(2*pi)
			expect_that(res, equals(c(3,2,1,8,5)))
			
			res <- frequenciesValidator(freq, N=10, steps=1:4)*10/(2*pi)
			expect_that(res, equals(c(3,2,1,2,5)))
			
			res <- frequenciesValidator(freq, N=10, steps=1:5)*10/(2*pi)
			expect_that(res, equals(c(3,2,1,5)))
			
			res <- frequenciesValidator(freq, N=10, steps=1:6)*10/(2*pi)
			expect_that(res, equals(c(1,2,3,5)))
			
			freq <- freq + 0.1 # No Fourier freq. any more!
			res <- (frequenciesValidator(freq, N=10, steps=6)-0.1)*10/(2*pi)
			expect_that(res, equals(c(2,3,5,8,11)))
			
			expect_that(res <- frequenciesValidator(freq, N=10, steps=3)*10/(2*pi),
					gives_warning())
			expect_that(res, equals(c(3,2,11,8,5)))
		}
)


test_that("clippedFT works as expected (in presence of ties)",{
      
      set.seed(2581)
      
      Y1 <- rnorm(64)
      Y1[31:35] <- Y1[1:5]
      #Y2 <- rank(Y1)/64
      Fn <- ecdf(Y1)
      Y2 <- Fn(Y1)
      
      freq <- 2*pi*(0:63)/64
      levels <- c(0.25,0.5)
      
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
