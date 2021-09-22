library(coda)                               # detailed MCMC output
library(parallel)                           # multithreading
library(runjags)                            # JAGS in R

load("df.rda")                              # green hawker population data
                                            # (externally available)
DATA <- with(df, list(Y = cbind(exuviae, egglaying_females),
                      Z = as.numeric(as.factor(area)),
                      X = cbind(rep(1, nrow(totals)),
                                as.numeric(treatment == "checker"),
                                as.numeric(year == "2016"),
                                as.numeric(year == "2017"),
                                as.numeric(year == "2018"),
                                emers_frac,
                                pH,
                                redox_V,
                                O2_frac,
                                EC_mS_cm,
                                temp_C,
                                water_depth_m,
                                sludge_thickness_m)))
DATA$ONES      <- rep(1, nrow(df))
DATA$NOBS      <- nrow(DATA$Y)
DATA$NCOV      <- ncol(DATA$X)
DATA$NREF      <- max(DATA$Z)
DATA$PRIORMU   <- rep(0, DATA$NCOV)
DATA$PRIORISIG <- diag(rep(0.1, DATA$NCOV))

MODEL <- '
model {                                     # Observation-specific models
    C <- 10000
    
    for (i in 1:NOBS) {
        ONES[i] ~ dbern(L[i] / C)
        
        L[i] <- P * equals(Y[i, 1], 0) * equals(Y[i, 1], 0) +
                    Q * equals(Y[i, 2], 0) * dpois(Y[i, 1], LA1[i]) +
                    R * equals(Y[i, 1], 0) * dpois(Y[i, 2], LA2[i]) +
                    (1 - P - Q - R) * LSUB[i]
        
        LSUB[i] <- sum(pow(exp(1), 0:min(Y[i, 1], Y[i, 2]) * log(LA0)
                                    + (Y[i, 1] - 0:min(Y[i, 1], Y[i, 2])) * log(LA1[i])
                                    + (Y[i, 2] - 0:min(Y[i, 1], Y[i, 2])) * log(LA2[i])
                                    - logfact(0:min(Y[i, 1], Y[i, 2]))
                                    - logfact(Y[i, 1] - 0:min(Y[i, 1], Y[i, 2]))
                                    - logfact(Y[i, 2] - 0:min(Y[i, 1], Y[i, 2]))
                                    - (LA0 + LA1[i] + LA2[i])))

                                        # Linear predictors
        log(LA1[i]) <- inprod(BETA_1, X[i, ]) + B_M[Z[i]]
        log(LA2[i]) <- inprod(BETA_2, X[i, ]) + B_N[Z[i]]
    }
    log(LA0) <- BETA_0

    for (j in 1:NREF) {                 # Random effects prior distributions
        B[j, 1:2] ~ dmnorm(c(0, 0), ISIGMA2)
        B_M[j] <- B[j, 1]
        B_N[j] <- B[j, 2]
    }
    
    OM1      ~  dgamma(0.5, 1E-2)       # Random effect variance MGH-t prior
    OM2      ~  dgamma(0.5, 1E-2)
    OM[1, 1] <- OM1
    OM[1, 2] <- 0
    OM[2, 1] <- 0
    OM[2, 2] <- OM2
    ISIGMA2  ~  dwish(inverse(OM), 3)
    
    
    P <- PI[1]                          # Zero-inflation models
    Q <- PI[2]
    R <- PI[3]
                                        # Fixed effects prior distributions
    for (k in 1:NCOV) {
        BETA_1[k] ~ dnorm(0, PREC1[k])
        BETA_2[k] ~ dnorm(0, PREC2[k])
        PREC1[k] ~ dgamma(2, 1)
        PREC2[k] ~ dgamma(2, 1)
    }
    BETA_0 ~ dnorm(0, PREC0)
    PREC0 ~ dgamma(2, 1)
    
    PI ~ ddirch(c(1, 1, 1, 1))          # Zero-inflation prior distribution
}'

PARAMS <- c('BETA_1', 'BETA_2', 'BETA_0',
            'P', 'Q', 'R',
            'B_M', 'B_N', 'ISIGMA2')

INITIAL <- replicate(8, list(list(.RNG.name = 'base::Mersenne-Twister',
                                  .RNG.seed = sample.int(n = 100000, size = 1))))

SAMPLE <- run.jags(model = MODEL, data = DATA, inits = INITIAL, monitor = PARAMS,
                   n.chains = 8, adapt = 1000, burnin = 10000, thin = 2500, sample = 100,
                   method = 'parallel', modules = c('ZIBGeometric', 'glm'),
                   factories = 'bugs::MNormal sampler off')
