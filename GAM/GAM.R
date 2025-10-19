GAM <- function(mu.link = "log", sigma.link = "log"){
  mstats <- checklink("mu.link", "GAM", substitute(mu.link),
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "GAM", substitute(sigma.link),
                      c("inverse", "log", "identity", "own"))
  structure(
    list(family = c("GAM", "Mi-Gamma"),
         parameters = list(mu=TRUE, sigma=TRUE),
         nopar = 2,
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)),
         mu.linkfun = mstats$linkfun,
         sigma.linkfun = dstats$linkfun,
         mu.linkinv = mstats$linkinv,
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta,
         sigma.dr = dstats$mu.eta,
         
         # First derivatives
         dldm = function(y,sigma,mu) (y-mu)/((sigma^2)*(mu^2)),
         dldd = function(y,mu,sigma) (2/sigma^3)*((y/mu)-log(y)+log(mu)+log(sigma^2)-1+digamma(1/(sigma^2))),
         
         # Second derivatives
         d2ldm2 = function(y, mu, sigma) {
           dldm <- (y-mu)/((sigma^2)*(mu^2))
           d2ldm2 <- - dldm * dldm
           d2ldm2
         },
         
         d2ldd2  = function(y, mu, sigma) {
           dldd <- (2/sigma^3)*((y/mu)-log(y)+log(mu)+log(sigma^2)-1+digamma(1/(sigma^2)))
           d2ldd2 <- - dldd * dldd
           d2ldd2
         },

         d2ldmdd = function(y, mu, sigma) {
           dldm <- (y-mu)/((sigma^2)*(mu^2))
           dldd <- (2/sigma^3)*((y/mu)-log(y)+log(mu)+log(sigma^2)-1+digamma(1/(sigma^2)))
           d2ldmdd <- - dldm * dldd
           d2ldmdd
         },
         
         # d2ldm2 = function(mu,sigma) -1/((sigma^2)*(mu^2)), 
         # d2ldd2 = function(sigma) (4/sigma^4)-(4/sigma^6)*trigamma((1/sigma^2)), 
         # d2ldmdd = function(y)  rep(0,length(y)),
         
         G.dev.incr = function(y,mu,sigma,...) -2*dGAM(y,mu,sigma,log=TRUE),
         rqres = expression(rqres(pfun="pGAM", type="Continuous",y=y,mu=mu,sigma=sigma)),
         
         mu.initial = expression({mu <- (y+mean(y))/2}),
         sigma.initial = expression({sigma <- rep(1,length(y))}) ,
         
         mu.valid = function(mu) all(mu > 0) ,
         sigma.valid = function(sigma) all(sigma > 0),
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

