# ZIBG-GLMM
Software package for ZIBG-GLMM Bayesian inference

## Quickstart
Install the JAGS `module` and run `example.R` to generate example data and generate a posterior parameter sample.

&nbsp;

## Contents
The contents of this package can be divided into three groups:

### Data generation
ZIBG distribution PMF  
`dzibg.cpp`  

MCMC algorithm to generate ZIBG-distributed data  
`mcmcsample.R`  

Functions defined in mcmcsample.R  
`rZIBG.rda`  

Script to generate data  
`generate.R`  

Data generated using `generate.R`  
`generated.rda`

&nbsp;


### JAGS model scripts invoked using `runjags`
Generate posterior sample for data in `generated.rda`  
`constantmodel.R`

Generate posterior sample for Aeshna viridis population data (externally available)  
`populationmodel.R`

Same as above, using the BZIP model in [1] instead  
`populationmodel_BZIP.R`

&nbsp;


### Custom JAGS module for ZIBG likelihood computations
Directory containing installation files
`module`

&nbsp;
&nbsp;



## References

[1] Majumdar, A., & Gries, C. (2010). Bivariate zero-inflated regression for count data: A Bayesian approach with application to plant counts. _The international journal of biostatistics, 6_(1).


