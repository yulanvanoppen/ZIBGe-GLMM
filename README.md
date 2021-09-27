# ZIBGe-GLMM
_Software package for Zero-Inflated Bivariate Geometric Generalized Linear Mixed Model (ZIBGe-GLMM) Bayesian inference_

&nbsp;



## Quickstart
Make sure the following software is installed and up-to-date:  
 - `R`  
 - `JAGS`  
 - `Boost`  

Also, make sure that the following `R` packages are installed and up-to-date:  
 - `runjags`  
 - `rjags`  
 - `Rcpp`  
 - `BH`  

Install the custom `/module` for `JAGS` (see `/module/README`) and run `/example.R` to generate example data and generate a posterior parameter sample.

&nbsp;



## Contents
The contents of this package can be divided into three groups:

### Data generation
ZIBGe distribution PMF  
`/dZIBGe.cpp`  

MCMC algorithm to generate ZIBGe-distributed data  
`/mcmcsample.R`  

Functions defined in `/mcmcsample.R`  
`/rZIBGe.rda`  

Script to generate data  
`/generate.R`  

Data generated using `generate.R`  
`/generated.rda`


### JAGS model scripts invoked using `runjags`
Generate posterior sample for data in `generated.rda`  
`/constantmodel.R`

Generate posterior sample for Aeshna viridis population data
`/populationmodel.R`

Same as above, using the BZIP model in [2] instead  
`/populationmodel_BZIP.R`  

Aeshna viridis population data loaded in `/populationmodel.R` and `/populationmodel_BZIP.R`  
`/populationdata.rda`


### Custom JAGS module for ZIBG likelihood computations
Directory containing installation files (see `module/README`)  
`/module`

&nbsp;



## Technical details

_Probability Mass Function_ (_PMF_) evaluations of the ZIBGe distribution need to be carefully implemented to avoid round-off errors when working with large multinomial coefficients. This problem is particularly pronounced when both components of the evaluated point are large. Therefore, multiple precision floating points are needed to store intermediate results since a `double` data type only supports precision up to 15 decimal digits. The PMF is implemented in a `C++` function (to be interfaced with `Rcpp` in `R`) to utilize the `cpp_bin_float` class from Boost's Multiprecision library (see [1]). Using high precision is computationally demanding, so to counteract this, powers already computed in each previous term are re-used to avoid evaluating the summands in the PMF directly. 

Likelihoods that cannot be computed using JAGS's built-in distributions are often dealt with using the zeros or ones trick (see [3], § 9.4). However, doing so prevents the use of multiple precision floating points. The modular character of `JAGS` makes it easy to extend the build-in distributions using custom (multivariate) distributions or sampling algorithms (see [4]). Analogous to `/dZIBGe.cpp`, the custom `/module` facilitates likelihood computations using high-precision intermediate computations.

Straightforward sampling methods, such as inverse transform sampling or rejection sampling, are inapplicable to the ZIBGe distribution; the former due to the nested sums appearing in its cumulative distribution function, and the latter because of the absence of a suitable proposal distribution. Therefore, a simple Metropolis-Hastings MCMC algorithm is used to generate ZIBG samples instead (see `/mcmcsample.R`). 

&nbsp;



## References

[1] Maddock, J., & Kormanyos, C. (2018). Boost multiprecision.

[2] Majumdar, A., & Gries, C. (2010). Bivariate zero-inflated regression for count data: A Bayesian approach with application to plant counts. _The International Journal of Biostatistics, 6_(1).

[3] Lunn, D., Jackson, C., Best, N., Thomas, A., & Spiegelhalter, D. (2012). _The BUGS Book: A Practical Introduction to Bayesian Analysis_ (CRC, Boca Raton, FL).

[4] Wabersich, D., & Vandekerckhove, J. (2014). Extending JAGS: A tutorial on adding custom distributions to JAGS (with a diffusion model example). _Behavior Research Methods, 46_(1), 15-28.

&nbsp;



## DISCLAIMER
_All software in this repository is covered by this disclaimer:_

While every effort is made to deliver high quality products, no guarantee is made that the products are free from defects. The software is provided “as is," and you use the software at your own risk.

No warranties are made as to performance, merchantability, fitness for a particular purpose, or any
other warranties whether expressed or implied.

No oral or written communication from or information provided by the author shall create a warranty.

Under no circumstances shall the author be liable for direct, indirect, special, incidental, or consequential damages resulting from the use, misuse, or inability to use this software, even if the author has been advised of the possibility of such damages.
