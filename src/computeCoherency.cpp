// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// computeSdNaive.cpp: Compute quantile-coherencies for smoothed
//					   quantile perioodgram
//
// Copyright (C) 2015 Tobias Kley
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

// int idx (int n, int d1, int k1, int d2, int k2, int b);

int idx (int n, int d1, int k1, int d2, int k2, int b, int N, int D1, int K1, int D2, int K2, int B) {
	int i = n;
	i += N * d1;
	i += N * D1 * k1;
	i += N * D1 * K1 * d2;
	i += N * D1 * K1 * D2 * k2;
	i += N * D1 * K1 * D2 * K2 * b * 1;
	return i;
}

//################################################################################
//' Workhorse function for \code{\link{getCoherency-SmoothedPG}}.
//'
//' C++ implementation to increase performance.
//'
//' @name .computeCoherency
//'
//' @keywords internals
//'
//' @param V a 3-dimensional array of complex numbers; dimensions are
//'          \code{[N, K1, K2]}, where \code{N} frequencies are
//'          \eqn{\omega_j := 2\pi j/N} for \eqn{j=0,\ldots,N}.
//' @param W a vector of length \code{W} of length \code{N} used for smoothing.
//'
//' @return Returns an array with complex numbers
//'         \eqn{\sigma(\tau_1, \tau_2, \omega_j} as defined in
//'         Kley et. al (2014), p. 26.
//'
//' @references
//' Barunik, J. & Kley, T. (2015+).
//' Quantile-Based Coherence: A Measure of Dependence for Econometric Variables.
//' \url{http://arxiv.org/abs/1401.8104}.
//################################################################################
// [[Rcpp::export(".computeCoherency")]]



ComplexVector computeCoherency(ComplexVector V, NumericVector d1, NumericVector d2) {

	// Dimensions of V
	IntegerVector D = V.attr("dim");

	if (D.size() == 4) {
		::Rf_error( "Cannot happen" ); // only good for multivariate time series!
	}
	int K1 = D(2);
	int K2 = D(4);
	int N = D(0);
	int B = D(5);

	int D1 = d1.size();
	int D2 = d2.size();

	// Declare a vector to hold the result
	ComplexVector res(N*D1*K1*D2*K2*B);
	res.attr("dim") = V.attr("dim");

	//Rcout << "Prod: " <<  N * D1 * K1 * D2 * K2 * B << std::endl;
	//Rcout << "B: " <<  B << std::endl;

	//try {

		for (int n = 0; n < N; n++) {
			for (int k1 = 0; k1 < K1; k1++) {
				for (int k2 = 0; k2 < K2; k2++) {
					for (int j1 = 0; j1 < D1; j1++) {
						for (int j2 = 0; j2 < D2; j2++) {
							for (int b = 0; b < B; b++) {

								//Rcout << "idx: " <<  idx(n, d1(j1)-1, k1, d2(j2)-1, k2, b, N, D(1), K1, D(3), K2, B) << " " << n << " " << d1(j1) << " " << k1 << " " << d2(j2) << " " << k2 << " " << b << " " << N << " " << D(1) << " " << K1 << " " << D(3) << " " << K2 << " " << B << std::endl;
								double nomR = V(idx(n, d1(j1)-1, k1, d2(j2)-1, k2, b, N, D(1), K1, D(3), K2, B)).r;
								double nomI = V(idx(n, d1(j1)-1, k1, d2(j2)-1, k2, b, N, D(1), K1, D(3), K2, B)).i;
								double den1 = V(idx(n, d1(j1)-1, k1, d1(j1)-1, k1, b, N, D(1), K1, D(3), K2, B)).r;
								double den2 = V(idx(n, d2(j2)-1, k2, d2(j2)-1, k2, b, N, D(1), K1, D(3), K2, B)).r;
								double fact = 1 / sqrt( den1 * den2 );
								res(idx(n, d1(j1)-1, k1, d2(j2)-1, k2, b, N, D(1), K1, D(3), K2, B)).r = nomR * fact;
								res(idx(n, d1(j1)-1, k1, d2(j2)-1, k2, b, N, D(1), K1, D(3), K2, B)).i = nomI * fact;
							}
						}
					}
				}
			}
		}

	//} catch( std::exception &ex ) {		// or use END_RCPP macro
	//	forward_exception_to_r( ex );
	//} catch(...) {
	//	::Rf_error( "c++ exception (unknown reason)" );
	//}
	return res;
}
