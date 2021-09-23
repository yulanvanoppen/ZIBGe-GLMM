library(coda)                           # detailed MCMC output
library(parallel)                       # multithreading
library(runjags)                        # invoke JAGS in R

load("df.rda")                          # Aeshna viridis population data
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
model {                                 # observation-specific models
    for (i in 1:NOBS) {
        Y[i, 1:2] ~ dzibg(MU[i], NU[i], TH, P, Q, R)

        MU[i] <- pow(pow(2 * (1 - P - R), 1 / (M[i] + 1)) - 1, -1)
        NU[i] <- pow(pow(2 * (1 - P - Q), 1 / (N[i] + 1)) - 1, -1)
        
        log(M[i])     <- inprod(BETA_M, X[i, ]) + B_M[Z[i]]
        log(N[i])     <- inprod(BETA_N, X[i, ]) + B_N[Z[i]]
    }
    logit(TH) <- BETA_TH                # correlation model

    for (i in 1:NREF) {                 # random effects prior distributions
        B[i, 1:2] ~ dmnorm(c(0, 0), ISIGMA2)
        B_M[i] <- B[i, 1]
        B_N[i] <- B[i, 2]
    }
    
    OM1      ~  dgamma(0.5, 1E-2)       # random effect variance MGH-t prior
    OM2      ~  dgamma(0.5, 1E-2)
    OM[1, 1] <- OM1
    OM[1, 2] <- 0
    OM[2, 1] <- 0
    OM[2, 2] <- OM2
    ISIGMA2  ~  dwish(inverse(OM), 3)
    
                                        # zero-inflation models
    P <- exp(BETA_P) / (1 + exp(BETA_P) + exp(BETA_Q) + exp(BETA_R)) / 2
    Q <- exp(BETA_Q) / (1 + exp(BETA_P) + exp(BETA_Q) + exp(BETA_R)) / 2
    R <- exp(BETA_R) / (1 + exp(BETA_P) + exp(BETA_Q) + exp(BETA_R)) / 2

                                        # fixed effects prior distributions
    BETA_M  ~ dmnorm(PRIORMU,       PRIORISIG)
    BETA_N  ~ dmnorm(PRIORMU,       PRIORISIG)
    
    BETA_TH ~ dnorm(0, 0.25)            # correlation prior distribution
    
    BETA_P  ~ dnorm(0, 0.25)            # zero-inflation prior distributions
    BETA_Q  ~ dnorm(0, 0.25)
    BETA_R  ~ dnorm(0, 0.25)
}'


PARAMS <- c('BETA_M', 'BETA_N', 'BETA_TH',
            'BETA_P', 'BETA_Q', 'BETA_R',
            'B_M', 'B_N', 'ISIGMA2')

INITIAL <- replicate(8, list(list(.RNG.name = 'base::Mersenne-Twister',
                                  .RNG.seed = sample.int(n = 100000, size = 1))))

SAMPLE <- run.jags(model = MODEL, data = DATA, inits = INITIAL, monitor = PARAMS,
                   n.chains = 8, adapt = 1000, burnin = 10000, thin = 2500, sample = 100,
                   method = 'parallel', modules = c('ZIBGeometric', 'glm'),
                   factories = 'bugs::MNormal sampler off')
