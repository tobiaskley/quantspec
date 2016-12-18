// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// computeSdNaive.cpp: Compute standard deviation for smoothed
//					   quantile perioodgram
//
// Copyright (C) 2014 Tobias Kley
//
// This file is part of quantspec.
//
// quantspec is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// quantspec is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include <cstdlib>
using namespace Rcpp;

//################################################################################
//' Workhorse function for \code{\link{initialize-ClippedFT}}.
//'
//' C++ implementation to increase performance.
//'
//' @name .generateIndMatrix
//'
//' @keywords internals
//'
//' @param Y time series to analysis; a T x D matrix.
//' @param pos.boot block bootstrap positions; a T x (B+1) matrix.
//' @param mult.boot multipliers for bootstrap; a T x (B+1) matrix.
//' @param levels levels; a vector of length K.
//' @param isRankBased passed from ClippedFT; boolean.
//'
//' @return Returns a matrix of dimension T x D*K*(B+1).
//################################################################################
// [[Rcpp::export(".generateIndMatrix")]]

NumericMatrix generateIndMatrix(NumericMatrix Y, NumericMatrix pos_boot, NumericMatrix mult_boot, NumericVector levels, bool isRankBased) {

	// Dimensions
	int T = Y.nrow();
	int D = Y.ncol();
	int B = pos_boot.ncol() - 1;
	int K = levels.size();
	

	// Declare a vector to hold the result
	NumericMatrix res( T, D*K*(B+1) );

	try {

   double FhatXt = 0;
   NumericVector Y_star(T);
   int idx = 0;
   
   for (int b = 0; b <= B; b++) {
      for (int d = 0; d < D; d++) {
         
         // determine block-bootstrapped sample
         for (int t = 0; t < T; t++) {
            Y_star[t] = Y(pos_boot(t, b), d);
         }
         
         for (int i = 0; i < K; i++) {
            for (int t = 0; t < T; t++) {
               
               // this seems something that should be done only once for every b/d?
               if (isRankBased) {
                  FhatXt = 0;
                  for (int j = 0; j < T; j++) {
                     if (Y_star[j] <= Y_star[t]) {
                        // FhatXt += mult_boot(j, b);
                        FhatXt += 1;
                     } 
                  }
                  FhatXt /= T;
               } else {
                  FhatXt = Y_star[t];
               } 
               
               idx = b*(K*D) + d*K + i;
               if (FhatXt <= levels[i]) {
                  // res(t, idx) = mult_boot(t, b);
                  res(t, idx) = mult_boot(t, b) * (1 - levels[i]);
               } else {
                  // res(t, idx) = 0;
                  res(t, idx) = mult_boot(t, b) * (0 - levels[i]);
               }
            }
         }
      }
   }

	} catch( std::exception &ex ) {		// or use END_RCPP macro
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "c++ exception (unknown reason)" );
	}
	return res;
}
