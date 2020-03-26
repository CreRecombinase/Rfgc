
// [[Rcpp::interfaces(r, cpp)]]


#include <Rcpp.h>
#include <fgc/fgc.h>
#include <RcppBlaze3.h>

#define RET_FOR_MEASURE(measure) \
    auto app = fgc::jsd::make_probdiv_applicator(X, measure); \
    auto ret = app.make_distance_matrix(measure, true); \
    return Rcpp::wrap(ret)

// [[Rcpp::export]]
Rcpp::NumericMatrix dist_matrixdd(blaze::DynamicMatrix<double> &X, int measure) {
    auto app = fgc::jsd::make_probdiv_applicator(X, measure);
    auto ret = app.make_distance_matrix(measure, true);
    return Rcpp::wrap(ret);
}
LLR  = 11
EMD  = 12
REVERSE_MKL  = 13
REVERSE_POISSON  = 14
UWLLR  = 15
TVD = TOTAL_VARIATION_DISTANCE  = 10
WASSERSTEIN=EMD  = 12
PSD = JSD  = 4
PSM = JSM  = 3

// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::LLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::LLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::LLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix llr_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::LLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix uwllr_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::UWLLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix uwllr_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::UWLLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix uwllr_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::UWLLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix uwllr_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::UWLLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(X, fgc::jsd::POISSON);
}
