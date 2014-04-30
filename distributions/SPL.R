# Shifted Power Law distribution for GAMLSS

SPL <- function(mu.link = "log", sigma.link = "logshiftto1")
{
  mstats <- checklink("mu.link", "Shifted Power Law", substitute(mu.link),
                      c("log", "identity", "own"))
  dstats <- checklink("sigma.link", "Shifted Power Law", substitute(sigma.link),
                      c("logshiftto1", "identity", "own"))

  structure(
    list(
      family        = c("SPL", "Shifted Power Law"),
      parameters    = list(mu = TRUE, sigma = TRUE),
      nopar         = 2,
      type          = "Continuous",
      mu.link       = as.character(substitute(mu.link)),
      sigma.link    = as.character(substitute(sigma.link)),
      mu.linkfun    = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
      mu.linkinv    = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
      mu.dr         = mstats$mu.eta,
      sigma.dr      = dstats$mu.eta,
      dldm          = function(y, mu, sigma) (sigma-1)/mu-sigma/(y+mu),
      d2ldm2        = function(y, mu, sigma) sigma/(y+mu)^2-(sigma-1)/mu^2,
      dldd          = function(y, mu, sigma) 1/(sigma-1)-log(y+mu)+log(mu),
      d2ldd2        = function(sigma) -1/(sigma-1)^2,
      d2ldmdd       = function(y, mu) 1/mu-1/(y+mu),
      G.dev.incr    = function(y, mu, sigma, ...) -2*dSPL(y, mu, sigma,
                                                          log = TRUE),
      rqres         = expression(rqres(pfun = "pSPL", type = "Continuous",
                                       y = y, mu = mu, sigma = sigma)),
      mu.initial    = expression(mu <- rep(5, length(y))),
      sigma.initial = expression(sigma <- rep(2.5, length(y))),
      mu.valid      = function(mu) all(mu > 0),
      sigma.valid   = function(sigma) all(sigma > 1),
      y.valid       = function(y) TRUE
    ),
    class = c("gamlss.family", "family")
  )
}

dSPL <- function(x, mu = 5, sigma = 2.5, log = FALSE)
{
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  d <- (sigma-1)/mu*((x+mu)/mu)^-sigma
  d[x < 0] <- 0

  if (log == TRUE) d <- log(d)

  d
}

pSPL <- function(q, mu = 5, sigma = 2.5, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  p <- 1-(mu/(q+mu))^(sigma-1)
  p[q < 0] <- 0

  if (lower.tail == FALSE) p <- 1-p
  if (log.p == TRUE) p <- log(p)

  p
}

qSPL <- function(p, mu = 5, sigma = 2.5, lower.tail = TRUE, log.p = FALSE)
{
  if (log.p == TRUE) p <- exp(p)
  if (lower.tail == FALSE) p <- 1-p

  if (any(p < 0) | any(p > 1)) stop("p must be between zero and one.")
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  q <- mu/(1-p)^(1/(sigma-1))-mu

  q
}

rSPL <- function(n, mu = 5, sigma = 2.5)
{
  if (any(n <= 0)) stop("n must be positive.")
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  n <- ceiling(n)
  p <- runif(n)
  r <- qSPL(p, mu = mu, sigma = sigma)

  r
}
