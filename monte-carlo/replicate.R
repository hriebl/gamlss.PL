replicate <- function(input, rep = 200) {
  n <- input$N

  if (min(input$y) >= 1) {
    gen <- "PL"
  } else {
    gen <- "SPL"
  }

  est <- input$family[1]

  c <- as.numeric(coef(input, "sigma")[1])
  beta1 <- as.numeric(coef(input, "sigma")[2])
  beta2 <- as.numeric(coef(input, "sigma")[3])

  for (i in 2:rep) {
    try(model <- simulate(n, gen, est))

    if (exists("model")) {
      c <- c(c, as.numeric(coef(model, "sigma")[1]))
      beta1 <- c(beta1, as.numeric(coef(model, "sigma")[2]))
      beta2 <- c(beta2, as.numeric(coef(model, "sigma")[3]))

      rm(model)
    } else {
      c <- c(c, NA)
      beta1 <- c(beta1, NA)
      beta2 <- c(beta2, NA)
    }
  }

  output <- list(n = n, data = data.frame(c, beta1, beta2))

  return(output)
}
