
// [[Rcpp::interfaces(r, cpp)]]


#include <Rcpp.h>
#include <fgc/fgc.h>
#include <RcppBlaze3.h>
using namespace fgc;

#define RET_FOR_MEASURE(measure) \
    auto app = fgc::jsd::make_probdiv_applicator(X, measure); \
    auto ret = app.make_distance_matrix(measure, true); \
    return Rcpp::wrap(ret)

// [[Rcpp::export]]
void display_constants() {
    fgc::jsd::detail::print_measures();
}
// [[Rcpp::export]]
void display_samplers() {
    std::fprintf(stderr, "0:\tfgc::coresets::BRAVERMAN_FELDMAN_LANG: BFL16\n"
                         "1:\tfgc::coresets::FELDMAN_LANGBERG: FL11\n"
                         "2:\tfgc::coresets::LUCIC_FAULKNER_KRAUSE_FELDMAN: LFKF17\n"
                         "3:\tfgc::coresets::VARADARAJAN_XIAO: VX11\n"
                         "4:\tfgc::coresets::LUCIC_FAULKNER_KRAUSE_FELDMAN: VX11\n");
    std::fprintf(stderr, "BFL is for generalizations of metrics with bicriteria approximations.\n"
                         "FL is for metric spaces and bicriteria approximations.\n"
                         "LFKF is for mixture modeling with bicriteria, while VX is for arbitrary metric spaces with constant-factor approximations\n"
                         "BKL (default) is for mu-similar divergences, as well as other metrics, using a bicriteria approximation. BKL is the default in Rfgc.\n");
}

template<typename MT>
Rcpp::List kmeans_coreset(MT &x, int k, size_t cs_size, uint64_t seed=0)
{
    std::mt19937_64 mt(seed);
    auto ret = fgc::coresets::kmeans_matrix_coreset(~x, k, mt, cs_size);
    return Rcpp::List::create(Rcpp::Named("weights") = Rcpp::wrap(std::move(ret.weights_)),
                              Rcpp::Named("points") = Rcpp::wrap(std::move(ret.mat_)),
                              Rcpp::Named("measure") = Rcpp::wrap("SQRL2"));
}

template<typename MT>
Rcpp::List construct_coreset(MT &x, int k, size_t cs_size, int measure, uint64_t seed=0, int sampler=(int)fgc::coresets::LBK) {
    fgc::jsd::ProbDivType type = (fgc::jsd::ProbDivType)measure;
    switch(type) {
        case fgc::jsd::TVD: case fgc::jsd::L1: std::fprintf(stderr, "Warning: D2 sampling is likely not sufficient for TVD and L1. Your coreset may have lower accuracy\n");
        default: ;
    }
    auto app = fgc::jsd::make_probdiv_applicator(~x, type, fgc::jsd::NONE);
    auto csampler = fgc::jsd::make_d2_coreset_sampler(app, k, seed, static_cast<blz::ElementType_t<MT> *>(nullptr), (fgc::coresets::SensitivityMethod)sampler);
    auto ics = csampler.sample(seed, cs_size);
    auto ret = fgc::coresets::index2matrix(ics, ~x);
    return Rcpp::List::create(Rcpp::Named("weights") = Rcpp::wrap(std::move(ret.weights_)),
                              Rcpp::Named("points") = Rcpp::wrap(std::move(ret.mat_)),
                              Rcpp::Named("measure") = Rcpp::wrap(fgc::jsd::detail::prob2str(type)));
}

// [[Rcpp::export]]
Rcpp::List construct_coresetdd(blaze::DynamicMatrix<double> &X, int k, size_t cs_size, int measure=1, int sampler=4, uint64_t seed=0) {
    return construct_coreset(X, k, cs_size, measure, seed);
}
// [[Rcpp::export]]
Rcpp::List construct_coresetdf(blaze::DynamicMatrix<float> &X, int k, size_t cs_size, int measure=1, int sampler=4, uint64_t seed=0) {
    return construct_coreset(X, k, cs_size, measure, seed);
}
// [[Rcpp::export]]
Rcpp::List construct_coresetsd(blaze::CompressedMatrix<double> &X, int k, size_t cs_size, int measure=1, int sampler=4, uint64_t seed=0) {
    return construct_coreset(X, k, cs_size, measure, seed);
}
// [[Rcpp::export]]
Rcpp::List construct_coresetsf(blaze::CompressedMatrix<float> &X, int k, size_t cs_size, int measure=1, int sampler=4, uint64_t seed=0) {
    return construct_coreset(X, k, cs_size, measure, seed);
}
// [[Rcpp::export]]
Rcpp::List kmeans_coresetdd(const blaze::DynamicMatrix<double> &X, int k, size_t cs_size, uint64_t seed=0) {
    return kmeans_coreset(X, k, cs_size, seed);
}
// [[Rcpp::export]]
Rcpp::List kmeans_coresetdf(const blaze::DynamicMatrix<float> &X, int k, size_t cs_size, uint64_t seed=0) {
    return kmeans_coreset(X, k, cs_size, seed);
}
// [[Rcpp::export]]
Rcpp::List kmeans_coresetsd(const blaze::CompressedMatrix<double> &X, int k, size_t cs_size, uint64_t seed=0) {
    return kmeans_coreset(X, k, cs_size, seed);
}
// [[Rcpp::export]]
Rcpp::List kmeans_coresetsf(const blaze::CompressedMatrix<float> &X, int k, size_t cs_size, uint64_t seed=0) {
    return kmeans_coreset(X, k, cs_size, seed);
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
