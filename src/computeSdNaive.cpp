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
//' Workhorse function for \code{\link{getSdNaive-SmoothedPG}}.
//'
//' C++ implementation to increase performance.
//'
//' @name .computeSdNaive
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
//' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2014).
//' Quantile Spectral Processes: Asymptotic Analysis and Inference.
//' \url{http://arxiv.org/abs/1401.8104}.
//################################################################################
// [[Rcpp::export(".computeSdNaive")]]

ComplexVector computeSdNaive(ComplexVector V, NumericVector W) {

	// Dimensions of V
	IntegerVector D = V.attr("dim");
	int K1 = D(1);
	int K2 = D(2);
	int N = W.size();

	// Declare a vector to hold the result
	ComplexVector res(Dimension(N,K1,K2));

	try {

	double S = 0;
	double *WW;
	WW = new double[N];
	for (int i = 0; i < N; i++) {
		WW[i] = W(i);
		S += WW[i];
	}

	double ***V1;
	double ***resRe;
	double ***resIm;
	double **V2;

	// Allocate memory for p
	V1 = new double**[K1];
	resRe = new double**[K1];
	resIm = new double**[K1];
	V2 = new double*[K1];
	for (int k1 = 0; k1 < K1; ++k1) {
		V2[k1] = new double[N];
		for (int s = 0; s < N; ++s) {
			V2[k1][s] = V(s+k1*D(0)+k1*D(0)*D(1)).r;
		}
		V1[k1] = new double*[K2];
		resRe[k1] = new double*[K2];
		resIm[k1] = new double*[K2];
		for (int k2 = k1; k2 < K2; ++k2) {
			V1[k1][k2-k1] = new double[N];
			resRe[k1][k2-k1] = new double[N];
			resIm[k1][k2-k1] = new double[N];
			for (int s = 0; s < N; ++s) {
				V1[k1][k2-k1][s] = (pow(V(s+k1*D(0)+k2*D(0)*D(1)).r,2) + pow(V(s+k1*D(0)+k2*D(0)*D(1)).i,2));
			}
		}
	}
	
	for (int j = 0; j < N; j++) {
		for (int k1 = 0; k1 < K1; k1++) {
			for (int k2 = k1; k2 < K2; k2++) {
			double z1 = 0;
			double z2 = 0;
			int jj = (N+j) % N;
			for (int s = 1; s < N; s++) {
				// first z1: for omega_j, z2: for -omega_j

				z1 += WW[(2*N+j-s) % N] * WW[(2*N+j-s) % N] * V2[k1][s] * V2[k2-k1][s];
				z1 += WW[(2*N+j-s) % N] * WW[(2*N+j+s) % N] * V1[k1][k2-k1][s];
				z2 += WW[(2*N+j-s) % N] * WW[(2*N-j-s) % N] * V2[k1][s] * V2[k2-k1][s];
				z2 += WW[(2*N+j-s) % N] * WW[(2*N-j+s) % N] * V1[k1][k2-k1][s];
				}
				if (k1 == k2) {
					resRe[k1][k2-k1][jj] = sqrt(std::max(0.0,z1)) / (S-WW[j % N]);
					resIm[k1][k2-k1][jj] = 0;
				} else {
					resRe[k1][k2-k1][jj] = sqrt(std::max(0.0,(z1+z2)/2)) / (S-WW[j % N]);
					resIm[k1][k2-k1][jj] = sqrt(std::max(0.0,(z1-z2)/2)) / (S-WW[j % N]);
				}
			}
		}
	}

	for (int k1 = 0; k1 < K1; ++k1) {
		for (int k2 = 0; k2 < 1; ++k2) {
			if (k1 == k2) {
				for (int s = 0; s < N; ++s) {
					res(s+N*k1+N*K1*k2).r = resRe[k1][0][s];
					res(s+N*k1+N*K1*k2).i = 0;
				}				
			} else {
				for (int s = 0; s < N; ++s) {
					if (k2 > k1) {
						res(s+N*k1+N*K1*k2).r = resRe[k1][k2-k1][s];
						res(s+N*k1+N*K1*k2).i = resIm[k1][k2-k1][s];
					} else {
						res(s+N*k1+N*K1*k2).r = resRe[k2][k1-k2][s];
						res(s+N*k1+N*K1*k2).i = resIm[k2][k1-k2][s];
					}
				}
			}
		}
	}
	
	// De-Allocate memory of V1, V2, W to prevent memory leak
	for (int k1 = 0; k1 < K1; ++k1) {
		for (int k2 = k1; k2 < K2; ++k2) {
			delete [] V1[k1][k2-k1];
			delete [] resRe[k1][k2-k1];
			delete [] resIm[k1][k2-k1];
		}
		delete [] V1[k1];
		delete [] resRe[k1];
		delete [] resIm[k1];
		delete [] V2[k1];
	}
	delete [] V1;
	delete [] resRe;
	delete [] resIm;
	delete [] V2;
	delete [] WW;

	} catch( std::exception &ex ) {		// or use END_RCPP macro
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "c++ exception (unknown reason)" );
	}
	return res;
}
