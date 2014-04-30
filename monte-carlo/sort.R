sort <- function(input, type = "est", mean = FALSE) {
  names <- c("c", "beta1", "beta2")

  est <- c()
  par <- c()
  n <- c()

  for (i in names) {
    if (i == "beta1") {
      true <- 0.2
    } else {
      true <- 0.1
    }

    for (j in input) {
      rep <- nrow(j$data)
      tmp <- j$data[, i]

      if (type == "bias") {
        tmp <- tmp-true
      }

      if (type == "mse") {
        tmp <- (tmp-true)**2
      }

      if (mean == TRUE) {
        for (k in rep:1) {
          tmp[k] <- mean(tmp[k:1])
        }
      }

      par <- c(par, rep(i, rep))
      n <- c(n, rep(j$n, rep))
      est <- c(est, tmp)
    }
  }

  output <- data.frame(par, n, est)

  return(output)
}
