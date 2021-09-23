#include <Rcpp.h>                                       // connect with R

#include <cmath>                                        // isfinite
#include <limits>                                       // most negative double
#include <boost/math/special_functions/binomial.hpp>    // binomial coefficient
#include <boost/multiprecision/cpp_bin_float.hpp>       // high precision floats

namespace mp = boost::multiprecision;                   // abbreviate for clarity
typedef mp::number<mp::cpp_bin_float<100> > mp_float;

// [[Rcpp::depends(BH)]]                                // library dependencies
// [[Rcpp::export]]                                     // make available in R
double dZIBGe(Rcpp::IntegerVector x, Rcpp::NumericVector par, bool log = false) {
    
    double y = x[0], z = x[1];                          // denote observed point and parameters
    double mu = par[0], nu = par[1], th = par[2], p = par[3], q = par[4], r = par[5];

    double numj = (1 - th) * mu * nu;                   // numerator and denominator factors
    double numy = mu + numj;
    double numz = nu + numj;
    double den = 1 + mu + nu + numj;
                                                        // start with term (coef*factor) j=0
    mp_float coef = boost::math::binomial_coefficient<mp_float>(y + z, y);
    mp_float factor = mp::pow(mp_float(numy), y) * mp::pow(mp_float(numz), z)
                                               / mp::pow(mp_float(den), y + z + 1);

    mp_float pmf = coef * factor;                        // bivariate geometric pmf
    
    for (int j = 0; j < std::min(y, z); ++j) {          // update term efficiently and add
        coef *= mp_float(y - j) * mp_float(z - j) / mp_float(y + z - j) / mp_float(j + 1);
        factor *= -numj / numy / numz * den;
        pmf += coef * factor;
    }
                                                        // add zero-inflation
    pmf = (1 - p - q - r) * pmf + p * (y == 0 && z == 0)
                                + q * (z == 0) * mp::pow(mp_float(mu / (mu + 1)), y) / (mu + 1)
                                + r * (y == 0) * mp::pow(mp_float(nu / (nu + 1)), z) / (nu + 1);
                                

    if (log)                                            // return (log) pmf
        return isfinite(mp::log(pmf).convert_to<double>()) ? mp::log(pmf).convert_to<double>()
                                                           : std::numeric_limits<double>::lowest();
    else
        return pmf.convert_to<double>();
}