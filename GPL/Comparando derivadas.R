
dldm_manual <- function(x, mu, sigma) {
  # constant term involving gamma functions
  C <- gamma(x + sigma) / (gamma(sigma) * gamma(x + 1))
  
  # numerator of the derivative inside the log
  num <- 2 * mu + C * (sigma * mu^(sigma - 1) * (mu + 1)^(1 - sigma) +
                         (1 - sigma) * mu^sigma * (mu + 1)^(-sigma))
  
  # denominator of that fraction (the argument inside the log)
  denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
  
  # derivative of log(p) with respect to mu
  dldm <- -(x + 2) / (mu + 1) + num / denom
  
  return(dldm)
}

dldm_compu = function(y, mu, sigma) {
  dm   <- gamlss::numeric.deriv(dGPL(y, mu, sigma, log=TRUE),
                                theta="mu",
                                delta=0.00001)
  dldm <- as.vector(attr(dm, "gradient"))
  dldm
}

y <- 0:10
mu <- 1.5
sigma <- 4.3

dldm_manual(y, mu, sigma)
dldm_compu(y, mu, sigma)



dldd_manual <- function(x, mu, sigma) {
  # constant term C(sigma)
  C <- gamma(x + sigma) / (gamma(sigma) * gamma(x + 1))
  
  # numerator inside the log derivative
  num <- C * mu^sigma * (mu + 1)^(1 - sigma) * 
    (digamma(x + sigma) - digamma(sigma) + log(mu / (mu + 1)))
  
  # denominator (the argument of the log)
  denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
  
  # derivative of log(p) with respect to sigma
  dldd <- num / denom
  
  return(dldd)
}

dldd_compu = function(y, mu, sigma) {
  dd   <- gamlss::numeric.deriv(dGPL(y, mu, sigma, log=TRUE),
                                theta="sigma",
                                delta=0.00001)
  dldd <- as.vector(attr(dd, "gradient"))
  dldd
}

y <- 0:10
mu <- 1.5
sigma <- 4.3

dldd_manual(y, mu, sigma)
dldd_compu(y, mu, sigma)



