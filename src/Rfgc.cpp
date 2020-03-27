
// [[Rcpp::interfaces(r, cpp)]]


#include <Rcpp.h>
#include <fgc/fgc.h>
#include <RcppBlaze3.h>

#define RET_FOR_MEASURE(measure) \
    auto app = fgc::jsd::make_probdiv_applicator(X, measure); \
    auto ret = app.make_distance_matrix(measure, true); \
    return Rcpp::wrap(ret)

// [[Rcpp::export]]
void display_constants() {
    fgc::jsd::detail::print_measures();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix dist_matrixdd(blaze::DynamicMatrix<double> &X, int arg) {
    fgc::jsd::ProbDivType measure = (fgc::jsd::ProbDivType) arg;
    RET_FOR_MEASURE(measure);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_matrixdf(blaze::DynamicMatrix<float> &X, int arg) {
    fgc::jsd::ProbDivType measure = (fgc::jsd::ProbDivType) arg;
    RET_FOR_MEASURE(measure);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_matrixsd(blaze::CompressedMatrix<double> &X, int arg) {
    fgc::jsd::ProbDivType measure = (fgc::jsd::ProbDivType) arg;
    RET_FOR_MEASURE(measure);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_matrixsf(blaze::CompressedMatrix<float> &X, int arg) {
    fgc::jsd::ProbDivType measure = (fgc::jsd::ProbDivType) arg;
    RET_FOR_MEASURE(measure);
}

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
    RET_FOR_MEASURE(fgc::jsd::UWLLR);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsd_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSD);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix jsm_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::JSM);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l2_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::L2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix l1_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::L1);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix sqrl2_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::SQRL2);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix mkl_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::MKL);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix tvd_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::TOTAL_VARIATION_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::POISSON);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_metric_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_METRIC);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_metric_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_METRIC);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_metric_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_METRIC);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_metric_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_METRIC);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_distance_matrixdd( blaze::DynamicMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_distance_matrixdf( blaze::DynamicMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_distance_matrixsd( blaze::CompressedMatrix<double> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_DISTANCE);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix bhattacharyya_distance_matrixsf( blaze::CompressedMatrix<float> &X) {
    RET_FOR_MEASURE(fgc::jsd::BHATTACHARYYA_DISTANCE);
}
