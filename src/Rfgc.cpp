
// [[Rcpp::interfaces(r, cpp)]]


#include <Rcpp.h>
#include <fgc/fgc.h>
#include <RcppBlaze3.h>

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
    return Rcpp::wrap(ret);
}
