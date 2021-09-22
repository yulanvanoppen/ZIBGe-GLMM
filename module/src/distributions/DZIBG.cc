#include <config.h>             // system configuration file, created by Autoconf
                                // and defined in configure.ac
#include "DZIBG.h"              // header file, containing class prototype
#include <rng/RNG.h>            // provides random functions
#include <util/nainf.h>         // provides na and inf functions etc.

#include <cmath>                // library for standard math operations
#include <algorithm>            // generic algorithms
#include <limits>               // represent the most negative double

#include <boost/math/special_functions/binomial.hpp>    // binomial coefficient
#include <boost/multiprecision/cpp_bin_float.hpp>       // high precision floats

using std::vector;              // vector is used in code
using std::min;                 // min is used in code
using std::max;                 // max is used in code

                                // high-precision (100 decimal places) floats
namespace mp = boost::multiprecision;
typedef mp::cpp_bin_float_100 mp_float; 

namespace jags::ZIBGeometric {  // module namespace

                                // constructor taking as arguments (1) the distribution
                                //   node's name as used in BUGS code and (2) the number
                                //   of parameters for that node
    DZIBG::DZIBG() : VectorDist("dzibg", 6) {}


                                // evaluate the log probability mass
    double DZIBG::logDensity(double const *x, unsigned int length, PDFType type,
		                      vector<double const *> const &parameters,
		                      vector<unsigned int> const &lengths,
		                      double const *lbound, double const *ubound) const {

        double y = x[0], z = x[1];                      // denote observed point and parameters
        double mu = *parameters[0], nu = *parameters[1], th = *parameters[2],
               pi1 = *parameters[3], pi2 = *parameters[4], pi3 = *parameters[5];

        double numj = (1 - th) * mu * nu;               // numerator and denominator factors
        double numy = mu + numj;
        double numz = nu + numj;
        double den = 1 + mu + nu + numj;

        mp_float coef = boost::math::binomial_coefficient<mp_float>(y + z, y);
        mp_float factor = mp::pow(mp_float(numy), y) * mp::pow(mp_float(numz), z)
                                                / mp::pow(mp_float(den), y + z + 1);

        mp_float pmf = coef * factor;                   // bivariate geometric pmf
        
        for (int j = 0; j < std::min(y, z); ++j) {      // update term efficiently and add
            coef *= mp_float(y - j) * mp_float(z - j) / mp_float(y + z - j) / mp_float(j + 1);
            factor *= -numj / numy / numz * den;
            pmf += coef * factor;
        }

        pmf = (1 - pi1 - pi2 - pi3) * pmf + pi1 * (y == 0 && z == 0)
                + pi2 * (z == 0) * mp::pow(mp_float(mu / (mu + 1)), y) / (mu + 1)
                + pi3 * (y == 0) * mp::pow(mp_float(nu / (nu + 1)), z) / (nu + 1);

        return jags_finite(mp::log(pmf).convert_to<double>()) ? mp::log(pmf).convert_to<double>()
                                                              : std::numeric_limits<double>::lowest();
    }

                                // draw random value - currently always (0, 0)
    void DZIBG::randomSample(double *x, unsigned int length,
			                    vector<double const *> const &parameters,
			                    vector<unsigned int> const &lengths, 
			                    double const *lbound, double const *ubound,
			                    RNG *rng) const {
        std::fill_n(x, 2, 0);
    }

                                // set equal to the mode - always (0, 0)
    void DZIBG::typicalValue(double *x, unsigned int length,
			                 vector<double const *> const &parameters,
			                 vector<unsigned int> const &lengths,
			                 double const *lbound, double const *ubound) const {
        std::fill_n(x, 2, 0);
    }

                                // set support for the distribution - (0, inf) x (0, inf)
    void DZIBG::support(double *lower, double *upper, unsigned int length,
			            vector<double const *> const &params, 
			            vector<unsigned int> const &lengths) const {
        std::fill_n(lower, 2, 0);
        std::fill_n(upper, 2, JAGS_POSINF);
    }

                                // returns true
    bool DZIBG::isSupportFixed(vector<bool> const &fixmask) const {
        return true;
    }

                                // verify that all parameters are scalar
    bool DZIBG::checkParameterLength (vector<unsigned int> const &parameters) const {
        return std::all_of(parameters.begin(),
                           parameters.end(),
                           [](unsigned int len) {return len == 1;});
    }

                                // verify positive means, correlation between -1 and 1, and
                                // probabilities that sum to 1
    bool DZIBG::checkParameterValue(vector<double const *> const &parameters,
			                        vector<unsigned int> const &lengths) const {
                                                        // denote observed point and parameters
        double mu = *parameters[0], nu = *parameters[1], th = *parameters[2],
                p = *parameters[3],  q = *parameters[4],  r = *parameters[5];

        return (mu >= 0 && nu >= 0 && th >= -1 && th <= 1 && p >= 0 &&  p <= 1 && 
                 q >= 0 &&  q <= 1 &&  r >=  0 &&  r <= 1 && p + q + r <= 1);
    }

                                // values are two-dimensional
    unsigned int DZIBG::length(vector<unsigned int> const &par) const {
        return 2;
    }
}