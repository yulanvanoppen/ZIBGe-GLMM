library(Rcpp)
sourceCpp('dzibg.cpp')                          # load ZIBG PMF

Q <- function(from, to) 2^-sum(!to)             # proposal distribution going Q(to | from)

propose <- function(current) {                  # sample neighbor from the proposal distribution
    move <- rnbinom(2, 1, mu=10) * sample(c(-1, 1), 2, replace=T)
    return(abs(current + move))
}

A <- function(from, to, pars, pmass)            # acceptance probability A(to | from)
    return( min(1, max(0, pmass(to, pars)) * Q(to, from) /
                          (pmass(from, pars) * Q(from, to))) )

mcmcsample <- function(n, pars, pmass, thin) { # generate sample using MCMC algorithm
    points <- matrix(0, n, 2)
    temp   <- matrix(0, thin, 2)

    for (point in 1:n) {                        # loop over n * thin iterations
        for (iter in 1:thin) {
            current <- temp[ifelse(iter == 1, thin, iter-1), ]
            proposal <- propose(current)        # propose new point
                                                
                                                # accept with probability A(proposal | current)
            if (A(current, proposal, pars, pmass) > runif(1))
                temp[iter, ] <- proposal
            else
                temp[iter, ] <- current
        }
        points[point, ] <- temp[thin, ]         # store every thin iterations
    }
    
    return(points)
}

rZIBG <- function(n, pars, thin=10000)          # generate ZIBG data
    mcmcsample(n, pars, pmass=dzibg, thin)

save(Q, propose, A, mcmcsample, rZIBG, file='rZIBG.rda')
