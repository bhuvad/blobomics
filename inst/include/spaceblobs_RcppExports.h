// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_spaceblobs_RCPPEXPORTS_H_GEN_
#define RCPP_spaceblobs_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace spaceblobs {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("spaceblobs", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("spaceblobs", "_spaceblobs_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in spaceblobs");
            }
        }
    }

    inline DataFrame fastHessian(NumericMatrix emat, NumericVector maskvec, IntegerVector x, IntegerVector y, int octaves, int threshold) {
        typedef SEXP(*Ptr_fastHessian)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_fastHessian p_fastHessian = NULL;
        if (p_fastHessian == NULL) {
            validateSignature("DataFrame(*fastHessian)(NumericMatrix,NumericVector,IntegerVector,IntegerVector,int,int)");
            p_fastHessian = (Ptr_fastHessian)R_GetCCallable("spaceblobs", "_spaceblobs_fastHessian");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_fastHessian(Shield<SEXP>(Rcpp::wrap(emat)), Shield<SEXP>(Rcpp::wrap(maskvec)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(octaves)), Shield<SEXP>(Rcpp::wrap(threshold)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<DataFrame >(rcpp_result_gen);
    }

    inline NumericMatrix smoothScales(DataFrame features, NumericVector maskvec, IntegerVector x, IntegerVector y) {
        typedef SEXP(*Ptr_smoothScales)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_smoothScales p_smoothScales = NULL;
        if (p_smoothScales == NULL) {
            validateSignature("NumericMatrix(*smoothScales)(DataFrame,NumericVector,IntegerVector,IntegerVector)");
            p_smoothScales = (Ptr_smoothScales)R_GetCCallable("spaceblobs", "_spaceblobs_smoothScales");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_smoothScales(Shield<SEXP>(Rcpp::wrap(features)), Shield<SEXP>(Rcpp::wrap(maskvec)), Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(y)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericMatrix >(rcpp_result_gen);
    }

}

#endif // RCPP_spaceblobs_RCPPEXPORTS_H_GEN_