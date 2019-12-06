// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// read_file_cpp
CharacterVector read_file_cpp(CharacterVector path);
RcppExport SEXP _rePROBE_read_file_cpp(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_file_cpp(path));
    return rcpp_result_gen;
END_RCPP
}
// read_file_cpp2
CharacterVector read_file_cpp2(std::string path);
RcppExport SEXP _rePROBE_read_file_cpp2(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_file_cpp2(path));
    return rcpp_result_gen;
END_RCPP
}
// reversecomp
CharacterVector reversecomp(CharacterVector x);
RcppExport SEXP _rePROBE_reversecomp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(reversecomp(x));
    return rcpp_result_gen;
END_RCPP
}
// prepare_vcf
void prepare_vcf(std::string file_in, std::string folder_out, std::string fs);
RcppExport SEXP _rePROBE_prepare_vcf(SEXP file_inSEXP, SEXP folder_outSEXP, SEXP fsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type file_in(file_inSEXP);
    Rcpp::traits::input_parameter< std::string >::type folder_out(folder_outSEXP);
    Rcpp::traits::input_parameter< std::string >::type fs(fsSEXP);
    prepare_vcf(file_in, folder_out, fs);
    return R_NilValue;
END_RCPP
}
// get_fasta
List get_fasta(CharacterVector path, List chrlist);
RcppExport SEXP _rePROBE_get_fasta(SEXP pathSEXP, SEXP chrlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type path(pathSEXP);
    Rcpp::traits::input_parameter< List >::type chrlist(chrlistSEXP);
    rcpp_result_gen = Rcpp::wrap(get_fasta(path, chrlist));
    return rcpp_result_gen;
END_RCPP
}
// read_fasta
CharacterVector read_fasta(CharacterVector path);
RcppExport SEXP _rePROBE_read_fasta(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_fasta(path));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rePROBE_read_file_cpp", (DL_FUNC) &_rePROBE_read_file_cpp, 1},
    {"_rePROBE_read_file_cpp2", (DL_FUNC) &_rePROBE_read_file_cpp2, 1},
    {"_rePROBE_reversecomp", (DL_FUNC) &_rePROBE_reversecomp, 1},
    {"_rePROBE_prepare_vcf", (DL_FUNC) &_rePROBE_prepare_vcf, 3},
    {"_rePROBE_get_fasta", (DL_FUNC) &_rePROBE_get_fasta, 2},
    {"_rePROBE_read_fasta", (DL_FUNC) &_rePROBE_read_fasta, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rePROBE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
