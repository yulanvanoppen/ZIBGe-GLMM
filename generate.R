library(Rcpp)
sourceCpp('dZIBGe.cpp')                         # load ZIBG PMF

load('rZIBGe.rda')                              # load ZIBG data-generating functions

M <- 4; N <- 2; th <- .5;                       # set distribution parameters
pi1 <- .1; pi2 <- .05; pi3 <- .025

mu <- ((2*(1 - pi1 - pi3))^(1/(M+1)) - 1)^-1    # reparameterize in terms of marginal means
nu <- ((2*(1 - pi1 - pi2))^(1/(N+1)) - 1)^-1

set.seed(0)
y <- rZIBGe(200, c(mu, nu, th, pi1, pi2, pi3))  # generate 200 data points

save(y, file='generated.rda')                   # store data