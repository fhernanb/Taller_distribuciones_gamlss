dGPL <- function(x, mu = 0.5, sigma = 0.5, log = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  # Ensure same length vector
  ly    <- max(length(x), length(mu), length(sigma))
  xx    <- rep(x, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  xx[x < 0] <- 0
  xx[is.infinite(x)] <- 1
  xx[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- 2 # No integers
  
  # pdf in log-scale
  term1 <- mu^2
  term2 <- (mu^sigma) * (mu + 1)^(1 - sigma) * gamma(x + sigma) / (gamma(sigma) * gamma(x + 1))
  p <- -(x + 2) * log(mu + 1) + log(term1 + term2)
  
  # Assign values for invalid x's
  p[x < 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  p[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- -Inf
  
  if (log == FALSE)
    p <- exp(p)
  return(p)
}
#' @export
#' @rdname dGPL
pGPL <- function(q, mu=0.5, sigma=1.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
  # Ensure same length vector
  ly    <- max(length(q), length(mu), length(sigma))
  qq    <- rep(q, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)
  
  # Temporal change for invalid x's
  qq[q < 0] <- 0
  qq[q == Inf] <- 0
  
  # For non-integer x's, the cumulative is the same as the lower integer
  qq <- as.integer(qq)
  
  # Auxiliary function
  fn <- function(q, mu, sigma) sum(dGPL(x=0:q, mu=mu, sigma=sigma))
  Vec_fn <- Vectorize(fn)
  
  # The cumulative
  cdf <- Vec_fn(q=qq, mu=mu, sigma=sigma)
  
  # Assign values for invalid x's
  cdf[q < 0] <- 0
  cdf[q == Inf] <- 1
  
  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)
  
  return(cdf)
}
#' @importFrom stats runif
#' @export
#' @rdname dGPL
rGPL <- function(n, mu = 0.5, sigma = 0.5) {
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  
  n <- ceiling(n)
  u <- runif(n=n)
  x <- qGPL(p=u, mu=mu, sigma=sigma)
  return(x)
}
#' @export
#' @rdname dGPL
qGPL <- function(p, mu=0.5, sigma=1.5, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  
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
  pp[p < 0]  <-  0.1
  pp[p > 1]  <-  0.1
  pp[p == 1] <-  0.1
  pp[p == 0] <-  0.1
  
  # Begin auxiliary function
  one_quantile_GPL <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      i <- 0
      prob <- dGPL(x=i, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      while (p >= F) {
        i <- i + 1
        prob <- dGPL(x=i, mu=mu, sigma=sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_GPL <- Vectorize(one_quantile_GPL)
  # End auxiliary function
  
  # The quantile
  q <- one_quantile_GPL(p=pp, mu=mu, sigma=sigma)
  
  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0
  
  return(q)
}





rGPL <- function(n, mu = 0.5, sigma = 0.5) {
  u <- runif(n, min=0, max=1)
  lambda <- ifelse(u < mu/(mu+1),
                   rexp(n=1, rate=mu),
                   rgamma(n=1, shape=sigma, rate=mu))
  
  y <- rpois(n=n, lambda=lambda)
  return(y)
}

rGPL <- function(n, mu = 0.5, sigma = 0.5) {
  u <- runif(n)
  lambda <- numeric(n)
  
  idx <- u < mu / (mu + 1)
  
  lambda[idx] <- rexp(n=sum(idx), rate = mu)
  lambda[!idx] <- rgamma(n=sum(!idx), shape = sigma, rate = mu)
  
  y <- rpois(n, lambda)
  return(y)
}

rGPL(n=2, mu=c(0.1, 0.2), sigma=c(2, 6))




