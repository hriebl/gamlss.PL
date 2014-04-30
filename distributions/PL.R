# Power Law distribution for GAMLSS

PL <- function(mu.link = "log", sigma.link = "logshiftto1")
{
  mstats <- checklink("mu.link", "Power Law", substitute(mu.link),
                      c("log", "identity", "own"))
  dstats <- checklink("sigma.link", "Power Law", substitute(sigma.link),
                      c("logshiftto1", "identity", "own"))

  structure(
    list(
      family        = c("PL", "Power Law"),
      parameters    = list(mu = FALSE, sigma = TRUE),
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
      dldm          = function(mu, sigma) (sigma-1)/mu,
      d2ldm2        = function(mu, sigma) (1-sigma)/mu^2,
      dldd          = function(y, mu, sigma) 1/(sigma-1)-log(y)+log(mu),
      d2ldd2        = function(sigma) -1/(sigma-1)^2,
      d2ldmdd       = function(mu) 1/mu,
      G.dev.incr    = function(y, mu, sigma, ...) -2*dPL(y, mu, sigma,
                                                         log = TRUE),
      rqres         = expression(rqres(pfun = "pPL", type = "Continuous",
                                       y = y, mu = mu, sigma = sigma)),
      mu.initial    = expression(mu <- rep(min(y), length(y))),
      sigma.initial = expression(sigma <- rep(2.5, length(y))),
      mu.valid      = function(mu) all(mu > 0),
      sigma.valid   = function(sigma) all(sigma > 1),
      y.valid       = function(y) TRUE
    ),
    class = c("gamlss.family", "family")
  )
}

dPL <- function(x, mu = 5, sigma = 2.5, log = FALSE)
{
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  d <- (sigma-1)/mu*(x/mu)^-sigma
  d[x < mu] <- 0

  if (log == TRUE) d <- log(d)

  d
}

pPL <- function(q, mu = 5, sigma = 2.5, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  p <- 1-(mu/q)^(sigma-1)
  p[q < mu] <- 0

  if (lower.tail == FALSE) p <- 1-p
  if (log.p == TRUE) p <- log(p)

  p
}

qPL <- function(p, mu = 5, sigma = 2.5, lower.tail = TRUE, log.p = FALSE)
{
  if (log.p == TRUE) p <- exp(p)
  if (lower.tail == FALSE) p <- 1-p

  if (any(p < 0) | any(p > 1)) stop("p must be between zero and one.")
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  q <- mu/(1-p)^(1/(sigma-1))

  q
}

rPL <- function(n, mu = 5, sigma = 2.5)
{
  if (any(n <= 0)) stop("n must be positive.")
  if (any(mu <= 0)) stop("mu must be positive.")
  if (any(sigma <= 1)) stop("sigma must be greater than one.")

  n <- ceiling(n)
  p <- runif(n)
  r <- qPL(p, mu = mu, sigma = sigma)

  r
}
