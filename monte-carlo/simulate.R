simulate <- function(n = 1000, gen = "PL", est = "PL", mu.start = NULL) {
  mu <- 1

  beta1 <- 0.2
  x1 <- rbinom(n, 1, 0.5)

  beta2 <- 0.1
  x2 <- runif(n, -2, 5)

  c <- 0.1

  eta <- beta1*x1+beta2*x2+c
  sigma <- exp(eta)+1

  if (gen == "PL") {
    y <- rPL(n, mu, sigma)
  } else {
    y <- rSPL(n, mu, sigma)
  }

  output <- gamlss(
    formula       = y~1,
    sigma.formula = y~x1+x2,
    family        = est,
    mu.start      = mu.start,
    control       = gamlss.control(gd.tol = Inf)
  )

  return(output)
}
