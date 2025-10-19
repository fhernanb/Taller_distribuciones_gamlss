# Example 1
# Generating some random values with
# known mu and sigma
set.seed(1234)
y <- rGPL(n=10000, mu=0.2, sigma=6)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=GPL,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ GPL
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(-1.6 + 5 * x1)
  sigma <- exp(1.7 - 5 * x2)
  y <- rGPL(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, mu=mu, sigma=sigma)
}

set.seed(12345)
datos <- gendat(n=10000)

mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=GPL, data=datos,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)

plot(mod2)
wp(mod2)
Rsq(mod2)


