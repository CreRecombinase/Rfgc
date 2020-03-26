// Copyright (C)  2017         Chingchuan Chen
// Copyright (C)  2010 - 2016  Dirk Eddelbuettel, Romain Francois and Douglas Bates
// Copyright (C)  2011         Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is based on RcppArmadillo.cpp and RcppEigen.h from RcppArmadillo and RcppEigen.
// This file is part of RcppBlaze3.
//
// RcppBlaze3.cpp: Rcpp/Blaze glue
//
// Copyright (C)  2017 Chingchuan Chen
//
// RcppBlaze3 is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppBlaze3 is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppBlaze3.  If not, see <http://www.gnu.org/licenses/>.


// [[Rcpp::interfaces(r, cpp)]


#include <RcppBlaze3.h>
#include <fgc/fgc.h>


// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixdd( blaze::DynamicMatrix<double> &X) {
    auto app = fgc::jsd::make_probdiv_applicator(X, fgc::jsd::LLR);
    auto ret = app.make_distance_matrix(fgc::jsd::LLR, true);
    return Rcpp::wrap(ret);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixdf( blaze::DynamicMatrix<float> &X) { 
    auto app = fgc::jsd::make_probdiv_applicator(X, fgc::jsd::LLR);
    auto ret = app.make_distance_matrix(fgc::jsd::LLR, true);
    return Rcpp::wrap(ret);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixsf( blaze::CompressedMatrix<float> &X) {
    auto app = fgc::jsd::make_probdiv_applicator(X, fgc::jsd::LLR);
    auto ret = app.make_distance_matrix(fgc::jsd::LLR, true);
    return Rcpp::wrap(ret);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixsd( blaze::CompressedMatrix<double> &X) {
    auto app = fgc::jsd::make_probdiv_applicator(X, fgc::jsd::LLR);
    auto ret = app.make_distance_matrix(fgc::jsd::LLR, true);
    return ret;
}
