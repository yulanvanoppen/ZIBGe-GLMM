library(Rcpp)
sourceCpp('dZIBGe.cpp')                          # load ZIBG PMF

load('rZIBGe.rda')                               # load ZIBG data-generating functions

set.seed(0)
y <- rZIBGe(200, c(5, 5, 0.5, 0.1, 0.05, 0.025)) # example: generate 200 data points

save(y, file='generated.rda')                   # store data