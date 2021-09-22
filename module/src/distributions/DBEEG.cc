#include <config.h>             // system configuration file, created by Autoconf
                                // and defined in configure.ac
#include "DBEEG.h"              // header file, containing class prototype
#include <rng/RNG.h>            // provides random functions
#include <util/nainf.h>         // provides na and inf functions etc.

#include <cmath>                // library for standard math operations
#include <algorithm>            // generic algorithms
#include <limits>               // represent the most negative double

#include <boost/math/special_functions/factorials.hpp>  // factorial
#include <boost/multiprecision/cpp_bin_float.hpp>       // high precision floats

using std::vector;              // vector is used in code

                                // high-precision (100 decimal places) floats
namespace mp = boost::multiprecision;
typedef mp::cpp_bin_float_100 myfloat;   

namespace jags::ZIBGeometric {  // module namespace

                                // constructor taking as arguments (1) the distribution
                                //   node's name as used in BUGS code and (2) the number
                                //   of parameters for that node
    DBEEG::DBEEG() : VectorDist("dzibeeg", 6) {}


                                // evaluate the log probability mass
    double DBEEG::logDensity(double const *x, unsigned int length, PDFType type,
		                      vector<double const *> const &parameters,
		                      vector<unsigned int> const &lengths,
		                      double const *lbound, double const *ubound) const {

        double y = x[0], z = x[1];  // rename for readability
        double th1 = *parameters[0], th2 = *parameters[1], b1 = *parameters[2],
                b2 = *parameters[3], rho = *parameters[4],  p = *parameters[5];

        double c1, c2, mu1, mu2, s1, s2, c11, c22;
        double pref1 = 1, pref2 = 1, pth1 = 1, pth2 = 1;
        for (std::size_t r = 1; r <= 10; ++r) {
            pref1 *= -(b1 - r + 1) / r;
            pref2 *= -(b2 - r + 1) / r;
            
            pth1 *= th1;
            pth2 *= th2;

            c1  += pref1 * (pth1 - 1) / (1 - pth1 * exp(-1));
            c2  += pref2 * (pth2 - 1) / (1 - pth2 * exp(-1));

            mu1 += -pref1 * pth1 / (1 - pth1);
            mu2 += -pref2 * pth2 / (1 - pth2);

            s1  += -pref1 * 2 * pow(pth1, 2) / pow(1 - pth1, 2);
            s2  += -pref2 * 2 * pow(pth2, 2) / pow(1 - pth2, 2);

            c11 += pref1 * pth1 * (pth1 - 1) / pow(1 - pth1 * exp(-1), 2);
            c22 += pref2 * pth2 * (pth2 - 1) / pow(1 - pth2 * exp(-1), 2);
        }

        s1 = sqrt(std::max(double(0), s1 - pow(mu1, 2)));
        s2 = sqrt(std::max(double(0), s2 - pow(mu2, 2)));

        double la = rho * s1 * s2 / (c11 - c1 * mu1) / (c22 - c2 * mu2);

        myfloat pmf = (1 + la * (exp(-y) - c1) * (exp(-z) - c2))
                          * (mp::pow(1 - mp::pow(myfloat(th1), y+1), b1)
                           - mp::pow(1 - mp::pow(myfloat(th1), y), b1))
                          * (mp::pow(1 - mp::pow(myfloat(th2), z+1), b2)
                           - mp::pow(1 - mp::pow(myfloat(th2), z), b2));

        pmf = (1 - p) * pmf + p * (y == 0 && z == 0);

        return jags_finite(mp::log(pmf).convert_to<double>()) ? mp::log(pmf).convert_to<double>()
                                                              : std::numeric_limits<double>::lowest();
    }

                                // draw random value - currently always (0, 0)
    void DBEEG::randomSample(double *x, unsigned int length,
			                    vector<double const *> const &parameters,
			                    vector<unsigned int> const &lengths, 
			                    double const *lbound, double const *ubound,
			                    RNG *rng) const {
        std::fill_n(x, 2, 0);
    }

                                // set equal to the mode - always (0, 0)
    void DBEEG::typicalValue(double *x, unsigned int length,
			                 vector<double const *> const &parameters,
			                 vector<unsigned int> const &lengths,
			                 double const *lbound, double const *ubound) const {
        std::fill_n(x, 2, 0);
    }

                                // set support for the distribution - (0, inf) x (0, inf)
    void DBEEG::support(double *lower, double *upper, unsigned int length,
			            vector<double const *> const &params, 
			            vector<unsigned int> const &lengths) const {
        std::fill_n(lower, 2, 0);
        std::fill_n(upper, 2, JAGS_POSINF);
    }

                                // returns true
    bool DBEEG::isSupportFixed(vector<bool> const &fixmask) const {
        return true;
    }

                                // verify that all parameters are scalar
    bool DBEEG::checkParameterLength (vector<unsigned int> const &parameters) const {
        return std::all_of(parameters.begin(),
                           parameters.end(),
                           [](unsigned int len) {return len == 1;});
    }

                                // verify positive means, correlation between -1 and 1, and
                                // probabilities that sum to 1
    bool DBEEG::checkParameterValue(vector<double const *> const &parameters,
			                        vector<unsigned int> const &lengths) const {
                                // rename for readability
        double th1 = *parameters[0], th2 = *parameters[1], b1 = *parameters[2],
                b2 = *parameters[3],  la = *parameters[4],  p = *parameters[5];

        return (th1 >= 0 && th1 <= 1 && th2 >= 0 && th2 <= 1 &&
                 b1 >= 0 &&  b2 >= 0 &&  la >= 0 &&   p >= 0 && p <= 1);
    }

                                // values are two-dimensional
    unsigned int DBEEG::length(vector<unsigned int> const &par) const {
        return 2;
    }



}