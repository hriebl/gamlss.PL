# Monte Carlo experiment to test the Power Law implementations for GAMLSS

library(gamlss)

source("../distributions/PL.R")
source("../distributions/SPL.R")
source("simulate.R")
source("replicate.R")
source("sort.R")

# Generate Power Law samples and estimate the corresponding GAMLSS. Plotting
# the biases of the estimates indicates that the results are consistent

set.seed(1337)
rep_250 <- replicate(simulate(250))

set.seed(1337)
rep_1000 <- replicate(simulate(1000))

set.seed(1337)
rep_5000 <- replicate(simulate(5000))

rep_all <- list(rep_250, rep_1000, rep_5000)
rep_sorted <- sort(rep_all, type = "bias")

par(mfrow = c(1, 3))

boxplot(
  rep_sorted$est[rep_sorted$n == 250 & rep_sorted$par == "c"],
  rep_sorted$est[rep_sorted$n == 1000 & rep_sorted$par == "c"],
  rep_sorted$est[rep_sorted$n == 5000 & rep_sorted$par == "c"],
  main = expression(beta[0])
)

axis(side = 1, at = c(1, 2, 3), labels = c(250, 1000, 5000))
abline(h = 0, col = "red", lwd = 2)

boxplot(
  rep_sorted$est[rep_sorted$n == 250 & rep_sorted$par == "beta1"],
  rep_sorted$est[rep_sorted$n == 1000 & rep_sorted$par == "beta1"],
  rep_sorted$est[rep_sorted$n == 5000 & rep_sorted$par == "beta1"],
  main = expression(beta[1])
)

axis(side = 1, at = c(1, 2, 3), labels = c(250, 1000, 5000))
abline(h = 0, col = "red", lwd = 2)

boxplot(
  rep_sorted$est[rep_sorted$n == 250 & rep_sorted$par == "beta2"],
  rep_sorted$est[rep_sorted$n == 1000 & rep_sorted$par == "beta2"],
  rep_sorted$est[rep_sorted$n == 5000 & rep_sorted$par == "beta2"],
  main = expression(beta[2])
)

axis(side = 1, at = c(1, 2, 3), labels = c(250, 1000, 5000))
abline(h = 0, col = "red", lwd = 2)

# Same procedure for the Shifted Power Law. The results are consistent again

set.seed(1337)
rep_250 <- replicate(simulate(250, "SPL", "SPL"), 50)

set.seed(1337)
rep_1000 <- replicate(simulate(1000, "SPL", "SPL"), 50)

set.seed(1337)
rep_5000 <- replicate(simulate(5000, "SPL", "SPL"), 50)

rep_all <- list(rep_250, rep_1000, rep_5000)
rep_sorted <- sort(rep_all, type = "bias")

boxplot(
  rep_sorted$est[rep_sorted$n == 250 & rep_sorted$par == "c"],
  rep_sorted$est[rep_sorted$n == 1000 & rep_sorted$par == "c"],
  rep_sorted$est[rep_sorted$n == 5000 & rep_sorted$par == "c"],
  main = expression(beta[0])
)

axis(side = 1, at = c(1, 2, 3), labels = c(250, 1000, 5000))
abline(h = 0, col = "red", lwd = 2)

boxplot(
  rep_sorted$est[rep_sorted$n == 250 & rep_sorted$par == "beta1"],
  rep_sorted$est[rep_sorted$n == 1000 & rep_sorted$par == "beta1"],
  rep_sorted$est[rep_sorted$n == 5000 & rep_sorted$par == "beta1"],
  main = expression(beta[1])
)

axis(side = 1, at = c(1, 2, 3), labels = c(250, 1000, 5000))
abline(h = 0, col = "red", lwd = 2)

boxplot(
  rep_sorted$est[rep_sorted$n == 250 & rep_sorted$par == "beta2"],
  rep_sorted$est[rep_sorted$n == 1000 & rep_sorted$par == "beta2"],
  rep_sorted$est[rep_sorted$n == 5000 & rep_sorted$par == "beta2"],
  main = expression(beta[2])
)

axis(side = 1, at = c(1, 2, 3), labels = c(250, 1000, 5000))
abline(h = 0, col = "red", lwd = 2)
