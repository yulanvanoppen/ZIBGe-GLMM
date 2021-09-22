library(coda)                           # detailed MCMC output
library(parallel)                       # multithreading
library(runjags)                        # JAGS in R

load('generated.rda')                   # (see /generate.R)
DATA <- list(Y = y, NOBS = nrow(y))

MODEL <- '
model {
    for (i in 1:NOBS) {                 # Observation-specific distribution
        Y[i, 1:2] ~ dzibg(MU, NU, TH, P, Q, R)
    }
                            # Mean reparameterization
    MU <- 1 / (pow(2 * (1 - P - R), 1 / (M + 1)) - 1)
    NU <- 1 / (pow(2 * (1 - P - Q), 1 / (N + 1)) - 1)
    
                                    # Marginal median and correlation models
    log(M)    <- BETA_M
    log(N)    <- BETA_N
    logit(TH) <- BETA_TH

                                        # Zero-inflation models
    P <- exp(BETA_P) / (1 + exp(BETA_P) + exp(BETA_Q) + exp(BETA_R)) / 2
    Q <- exp(BETA_Q) / (1 + exp(BETA_P) + exp(BETA_Q) + exp(BETA_R)) / 2
    R <- exp(BETA_R) / (1 + exp(BETA_P) + exp(BETA_Q) + exp(BETA_R)) / 2
    
                                        # Parameter prior distributions
    BETA_M  ~ dnorm(0, 0.1)
    BETA_N  ~ dnorm(0, 0.1)
    BETA_TH ~ dnorm(0, 0.25)
    
    BETA_P  ~ dnorm(0, 0.25)
    BETA_Q  ~ dnorm(0, 0.25)
    BETA_R  ~ dnorm(0, 0.25)
    
}'

PARAMS <- c('BETA_M', 'BETA_N', 'BETA_TH', 'BETA_P', 'BETA_Q', 'BETA_R')

INITIAL <- replicate(8, list(list(.RNG.name = 'base::Mersenne-Twister',
                                  .RNG.seed = sample.int(n = 100000, size = 1))))

SAMPLE <- run.jags(model = MODEL, inits = INITIAL, data = DATA, n.chains = 8,
                   monitor = PARAMS, adapt = 100, burnin = 50, thin = 5, sample = 100,
                   method = 'parallel', modules = c('ZIBGeometric', 'glm'))






