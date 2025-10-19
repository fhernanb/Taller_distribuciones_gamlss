dGAM <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  
  # Ensure same length vector
  ly    <- max(length(x), length(mu), length(sigma))
  xx    <- rep(x, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  xx[x <= 0] <- 0.5
  xx[is.infinite(x)] <- 0.5
  
  # pdf in log-scale
  p <- (1/sigma^2)*log(xx/(mu*sigma^2))-xx/(mu*sigma^2)-log(xx)-lgamma(1/sigma^2)
  
  # Assign values for invalid x's
  p[x <= 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  
  if (log == FALSE)
    p <- exp(p)
  
  return(p)
}

pGAM <- function(q, mu=1, sigma=1, lower.tail=TRUE, log.p=FALSE){
  
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  # Ensure same length vector
  ly    <- max(length(q), length(mu), length(sigma))
  qq    <- rep(q, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  qq[q <= 0] <- 0.5
  qq[q == Inf] <- 0.5
  
  # The cumulative
  p1 <- zipfR::Igamma(a=1/sigma^2, x=qq/(mu*sigma^2))
  p2 <- gamma(1/sigma^2)
  cdf <- p1 / p2
  
  # Assign values for invalid x's
  cdf[q <= 0] <- 0
  cdf[q == Inf] <- 1

  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)
  
  return(cdf)
}


qGAM <- function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", ""))
  
  # To adjust the probability
  if (log.p == TRUE)
    p <- exp(p)
  if (lower.tail == FALSE)
    p <- 1 - p
  
  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp <- rep(p, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid p's
  pp[p < 0]  <-  0.5
  pp[p > 1]  <-  0.5
  pp[p == 1] <-  0.5
  pp[p == 0] <-  0.5
  
  # The quantile
  q <- qgamma(pp, shape=1/sigma^2, scale=mu*sigma^2, 
              lower.tail=lower.tail)
  
  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0
  
  return(q)
}

rGAM <- function(n, mu=1, sigma=1){
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  
  n <- ceiling(n)
  u <- runif(n=n)
  x <- qGAM(p=u, mu=mu, sigma=sigma)
  return(x)
}


