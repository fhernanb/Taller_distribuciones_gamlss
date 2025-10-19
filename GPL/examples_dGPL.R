# Example 1
# Plotting the mass function for different parameter values

x_max <- 50
probs1 <- dGPL(x=0:x_max, mu=0.1, sigma=2)
probs2 <- dGPL(x=0:x_max, mu=0.5, sigma=5)
probs3 <- dGPL(x=0:x_max, mu=0.2, sigma=6)

# To plot the first k values
plot(x=0:x_max, y=probs1, type="o", lwd=2, col="green4", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for GPL",
     ylim=c(0, 0.040))
points(x=0:x_max, y=probs2, type="o", lwd=2, col="tomato")
points(x=0:x_max, y=probs3, type="o", lwd=2, col="black")
legend("topright", col=c("green4", "tomato", "black"), lwd=3,
       legend=c("mu=0.001, sigma=0.52 ",
                "mu=0.001, sigma=0.85",
                "mu=0.001, sigma=1.5"))

# Example 2
# Checking if the cumulative curves converge to 1

x_max <- 50
cumulative_probs1 <- pGPL(q=0:x_max, mu=0.1, sigma=2)
cumulative_probs2 <- pGPL(q=0:x_max, mu=0.5, sigma=5)
cumulative_probs3 <- pGPL(q=0:x_max, mu=0.6, sigma=6)

plot(x=0:x_max, y=cumulative_probs1, col="green4",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative probability for GPL",
     xlab="X", ylab="Probability")
points(x=0:x_max, y=cumulative_probs2, type="o", col="tomato")
points(x=0:x_max, y=cumulative_probs3, type="o", col="black")
legend("bottomright", col=c("green4", "tomato", "black"), lwd=3,
       legend=c("mu=0.001, sigma=0.52 ",
                "mu=0.001, sigma=0.85",
                "mu=0.001, sigma=1.5"))

# Example 3
# Comparing the random generator output with the theoretical probabilities

x_max <- 150
mu <- 0.1
sigma <- 2
probs1 <- dGPL(x=0:x_max, mu=mu, sigma=sigma)
names(probs1) <- 0:x_max

x <- rGPL(n=10000, mu=mu, sigma=sigma)
probs2 <- prop.table(table(x))

cn <- union(names(probs1), names(probs2))
height <- rbind(probs1[cn], probs2[cn])
nombres <- cn
mp <- barplot(height, beside = TRUE, names.arg = nombres,
              col=c("dodgerblue3","firebrick3"), las=1,
              xlab="X", ylab="Proportion")
legend("topright",
       legend=c("Theoretical", "Simulated"),
       bty="n", lwd=3,
       col=c("dodgerblue3","firebrick3"), lty=1)

# Example 4
# Checking the quantile function

mu <- 0.2
sigma <- 6
p <- seq(from=0, to=1, by=0.01)
qxx <- qGPL(p=p, mu=mu, sigma=sigma, lower.tail=TRUE, log.p=FALSE)
plot(p, qxx, type="s", lwd=2, col="green3", ylab="quantiles",
     main="Quantiles of GPL(mu = sigma = 0.2)")
