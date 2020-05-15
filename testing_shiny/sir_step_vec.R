sir_step <- function (sims, S, E, I1, I2, R, N, beta, theta,gamma_I1I2,gamma_I2R, delta.t, ...) {
  # adapted from Aaron King's code
  
  dN_SE <- rbinom(n=sims,size=S,prob=1-exp(-beta*(I1+I2)/N*delta.t))
  dN_EI1 <- rbinom(n=sims,size=E,prob=1-exp(-theta*delta.t))
  dN_I1I2 <- rbinom(n=sims,size=I1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R <- rbinom(n=sims,size=I2,prob=1-exp(-gamma_I2R*delta.t))
  #browser()
  S <- S - dN_SE
  E <- E + dN_SE - dN_EI1
  I1 <- I1 + dN_EI1 - dN_I1I2
  I2 <- I2 + dN_I1I2 - dN_I2R
  R <- R + dN_I2R
  cbind(S ,  E,  I1,  I2, R, dN_I1I2 ) # assume that I1->I2 is when cases become detectable
}
