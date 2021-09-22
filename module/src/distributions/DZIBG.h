#ifndef DZIBG_H_
#define DZIBG_H_
#include <distribution/VectorDist.h>    // JAGS vector distribution base class

namespace jags::ZIBGeometric {
                                // vector distribution class
    class DZIBG : public VectorDist {   
        public:
            DZIBG();            // constructor

                                // evaluate the log probability mass
            double logDensity(double const *x, unsigned int length, PDFType type,
		                      std::vector<double const *> const &parameters,
		                      std::vector<unsigned int> const &lengths,
		                      double const *lbound, double const *ubound) const;

                                // draw random value - currently always (0, 0)
            void randomSample(double *x, unsigned int length,
			                     std::vector<double const *> const &parameters,
			                     std::vector<unsigned int> const &lengths, 
			                     double const *lbound, double const *ubound,
			                     RNG *rng) const;
            
                                // set equal to the mode - always (0, 0)
            void typicalValue(double *x, unsigned int length,
			                  std::vector<double const *> const &parameters,
			                  std::vector<unsigned int> const &lengths,
			                  double const *lbound, double const *ubound) const;

                                // set support for the distribution - (0, inf) x (0, inf)
            void support(double *lower, double *upper, unsigned int length,
			             std::vector<double const *> const &params, 
			             std::vector<unsigned int> const &lengths) const;
            
                                // returns true
            bool isSupportFixed(std::vector<bool> const &fixmask) const;

                                // verify that parameters are scalar
            bool checkParameterLength(std::vector<unsigned int> const &parameters) const;

                                // verify positive means, correlation between -1 and 1, and
                                // probabilities that sum to 1
            bool checkParameterValue(std::vector<double const *> const & parameters,
                                     std::vector<unsigned int> const &lengths) const;

                                // values are two-dimensional
            unsigned int length (std::vector<unsigned int> const &par) const;
    };

}

#endif                          // DZIBG_H_