GPL <- function(mu.link="log", sigma.link="log") {
  
  mstats <- checklink("mu.link", "GPL",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "GPL",
                      substitute(sigma.link), c("log"))
  
  structure(list(family=c("GPL", "Poisson-Generalized Lindley"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Discrete",
                 
                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 
                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 
                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 
                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 
                 # First derivatives
                 
                 dldm = function(y, mu, sigma) {
                   C <- gamma(y + sigma) / (gamma(sigma) * gamma(y + 1))
                   num <- 2 * mu + C * (sigma * mu^(sigma - 1) * (mu + 1)^(1 - sigma) +
                                          (1 - sigma) * mu^sigma * (mu + 1)^(-sigma))
                   denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
                   dldm <- -(y + 2) / (mu + 1) + num / denom
                   
                   dldm
                 },
                 
                 dldd = function(y, mu, sigma) {
                   C <- gamma(y + sigma) / (gamma(sigma) * gamma(y + 1))
                   num <- C * mu^sigma * (mu + 1)^(1 - sigma) * 
                     (digamma(y + sigma) - digamma(sigma) + log(mu / (mu + 1)))
                   denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
                   dldd <- num / denom
                   
                   dldd
                 },
                 
                 # Second derivatives
                 
                 d2ldm2 = function(y, mu, sigma) {
                   C <- gamma(y + sigma) / (gamma(sigma) * gamma(y + 1))
                   num <- 2 * mu + C * (sigma * mu^(sigma - 1) * (mu + 1)^(1 - sigma) +
                                          (1 - sigma) * mu^sigma * (mu + 1)^(-sigma))
                   denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
                   dldm <- -(y + 2) / (mu + 1) + num / denom
                   
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },
                 
                 d2ldmdd = function(y, mu, sigma) {
                   C <- gamma(y + sigma) / (gamma(sigma) * gamma(y + 1))
                   num <- 2 * mu + C * (sigma * mu^(sigma - 1) * (mu + 1)^(1 - sigma) +
                                          (1 - sigma) * mu^sigma * (mu + 1)^(-sigma))
                   denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
                   dldm <- -(y + 2) / (mu + 1) + num / denom
                   
                   C <- gamma(y + sigma) / (gamma(sigma) * gamma(y + 1))
                   num <- C * mu^sigma * (mu + 1)^(1 - sigma) * 
                     (digamma(y + sigma) - digamma(sigma) + log(mu / (mu + 1)))
                   denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
                   dldd <- num / denom

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },
                 
                 d2ldd2  = function(y, mu, sigma) {
                   C <- gamma(y + sigma) / (gamma(sigma) * gamma(y + 1))
                   num <- C * mu^sigma * (mu + 1)^(1 - sigma) * 
                     (digamma(y + sigma) - digamma(sigma) + log(mu / (mu + 1)))
                   denom <- mu^2 + C * mu^sigma * (mu + 1)^(1 - sigma)
                   dldd <- num / denom
                   
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },
                 
                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dGPL(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pGPL", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),
                 
                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_GPL(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_GPL(y)[2], length(y)) ),
                 
                 #mu.initial    = expression(mu    <- rep(1, length(y)) ),
                 #sigma.initial = expression(sigma <- rep(1, length(y)) ),
                 
                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 
                 y.valid = function(y) all(y >= 0)
                 
  ),
  class=c("gamlss.family", "family"))
}
#'
#' Initial values for GPL
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_GPL <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_GPL,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' logLik function for GPL
#' @description Calculates logLik for GPL distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @return returns the loglikelihood given the parameters and random sample.
#' @keywords internal
#' @export
logLik_GPL <- function(logparam=c(0, 0), x){
  return(sum(dGPL(x,
                     mu    = exp(logparam[1]),
                     sigma = exp(logparam[2]),
                     log=TRUE)))
}