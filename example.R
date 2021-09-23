library(coda)                               # detailed MCMC output
library(parallel)                           # multithreading
library(Rcpp)                               # compile C++ files
library(runjags)                            # invoke JAGS in R

sourceCpp('dZIBGe.cpp')                      # load ZIBG PMF

load('rZIBGe.rda')                           # load ZIBG data-generating functions

set.seed(0)                                 # fix random number generator
y <- rZIBGe(200, c(5, 5, 0.5, 0.1, 0, 0))    # generate 200 data points
plot(y)

DATA <- list(Y = y, NOBS = nrow(y))         # prepare data for JAGS

MODEL <- '
model {
    for (i in 1:NOBS) {                     # observation-specific distribution
        Y[i, 1:2] ~ dzibg(MU, NU, TH, P, 0, 0)
    }
                                            # mean reparameterization
    MU <- 1 / (pow(2 * (1 - P), 1 / (M + 1)) - 1)
    NU <- 1 / (pow(2 * (1 - P), 1 / (N + 1)) - 1)
    
    M  ~ dlnorm(0, 0.1)                     # prior distributions
    N  ~ dlnorm(0, 0.1)
    TH ~ dunif(0, 1)
    P  ~ dunif(0, .5)
    
}'

PARAMS <- c('M', 'N', 'TH', 'P')           # monitored parameters

                                            # chain seeds
INITIAL <- replicate(8, list(list(.RNG.name = 'base::Mersenne-Twister',
                                  .RNG.seed = sample.int(n = 100000, size = 1))))

                                            # generate posterior sample
                                            # (make sure the 'ZIBGeometric' module is installed)
SAMPLE <- run.jags(model = MODEL, inits = INITIAL, data = DATA, n.chains = 8,
                   monitor = PARAMS, adapt = 100, burnin = 50, thin = 5, sample = 100,
                   method = 'parallel', modules = c('ZIBGeometric', 'glm'))

print(summary(SAMPLE))






