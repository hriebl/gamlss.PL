# Analysis of German city sizes with the Power Law distributions for GAMLSS

# Prepare data

cities <- read.csv("cities.csv")

states <- c("Berlin", "Bremen", "Hamburg")
cities <- cities[! cities$state %in% states, ]

# Plot histograms on a linear scale and on a log-log scale. The red line in the
# right plot seems a pretty good fit, thus 2000 could be the lower bound of the
# Power Law behavior

breaks <- (max(cities$total)-min(cities$total))/100
hist <- hist(cities$total, breaks, plot = FALSE)

par(mfrow = c(1, 2))

plot(
  x    = hist$mids[1:500],
  y    = hist$counts[1:500],
  type = "l",
  main = "Linear Scale",
  xlab = "City Size",
  ylab = "Count",
)

plot(
  x    = hist$mids[1:500],
  y    = hist$counts[1:500],
  type = "l",
  main = "Log-Log Scale",
  xlab = "City Size",
  ylab = "Count",
  log  = "xy"
)

x <- hist$mids[20:500][hist$counts[20:500] != 0]
y <- hist$counts[20:500][hist$counts[20:500] != 0]
coef <- lm(log10(y) ~ log10(x))$coefficients

abline(coef, col = "red", lwd = 2)

# Look for optimal lower bound minimizing the KS-statistic. The calculation
# takes a long time, so we commented it out. According to the KS-statistic,
# 17487 is the optimal lower bound

cities_ks <- read.csv("cities-ks.csv")

par(mfrow = c(1, 1))

plot(
  x    = cities_ks$xmin,
  y    = cities_ks$ks,
  type = "l",
  xlab = expression(x[min]),
  ylab = "KS-statistic"
)

# library(poweRlaw)

# m <- displ(cities$total)
# est <- estimate_xmin(m)

# Apply different GAMLSS distributions to the city sizes without a lower bound.
# The Pareto distribution and the Shifted Power Law fit best because they have
# the lowest AIC

library(gamlss)

source("../distributions/PL.R")
source("../distributions/SPL.R")

results <- as.data.frame(matrix(nrow = 0, ncol = 3))
names(results) <- c("GAMLSS Distribution", "Estimate for Farming", "AIC")

cities_no <- gamlss(
  formula       = total ~ farming,
  sigma.formula = total ~ 1,
  family        = NO,
  data          = cities
)

cities_pa <- gamlss(
  formula       = total ~ 1,
  sigma.formula = total ~ farming,
  data          = cities,
  family        = PARETO2o,
  sigma.start   = 2.5,
  control       = gamlss.control(n.cyc = Inf)
)

cities_pl <- gamlss(
  formula       = total ~ 1,
  sigma.formula = total ~ farming,
  family        = PL,
  data          = cities
)

cities_spl <- gamlss(
  formula       = total ~ 1,
  sigma.formula = total ~ farming,
  family        = SPL,
  data          = cities,
  control       = gamlss.control(n.cyc = Inf)
)

info <- c(round(coef(cities_no, "mu")[2], 3), round(AIC(cities_no), 3))
results[(nrow(results) +  1), ] <- c("Normal distribution", info)

info <- c(round(coef(cities_pa, "sigma")[2], 3), round(AIC(cities_pa), 3))
results[(nrow(results) +  1), ] <- c("Pareto distribution", info)

info <- c(round(coef(cities_pl, "sigma")[2], 3), round(AIC(cities_pl), 3))
results[(nrow(results) +  1), ] <- c("Power Law distribution", info)

info <- c(round(coef(cities_spl, "sigma")[2], 3), round(AIC(cities_spl), 3))
results[(nrow(results) +  1), ] <- c("Shifted Power Law distribution", info)

results

# Same procedure as before. With 2000 as the lower bound, the Power Law
# distribution turns out to have the lowest AIC

cities <- cities[cities$total > 2000, ]

results <- as.data.frame(matrix(nrow = 0, ncol = 3))
names(results) <- c("GAMLSS Distribution", "Estimate for Farming", "AIC")

cities_no <- gamlss(
  formula       = total ~ farming,
  sigma.formula = total ~ 1,
  family        = NO,
  data          = cities
)

cities_pa <- gamlss(
  formula       = total ~ 1,
  sigma.formula = total ~ farming,
  data          = cities,
  family        = PARETO2o,
  sigma.start   = 2.5,
  control       = gamlss.control(n.cyc = Inf)
)

cities_pl <- gamlss(
  formula       = total ~ 1,
  sigma.formula = total ~ farming,
  family        = PL,
  data          = cities
)

cities_spl <- gamlss(
  formula       = total ~ 1,
  sigma.formula = total ~ farming,
  family        = SPL,
  data          = cities,
  control       = gamlss.control(n.cyc = Inf)
)

info <- c(round(coef(cities_no, "mu")[2], 3), round(AIC(cities_no), 3))
results[(nrow(results) +  1), ] <- c("Normal distribution", info)

info <- c(round(coef(cities_pa, "sigma")[2], 3), round(AIC(cities_pa), 3))
results[(nrow(results) +  1), ] <- c("Pareto distribution", info)

info <- c(round(coef(cities_pl, "sigma")[2], 3), round(AIC(cities_pl), 3))
results[(nrow(results) +  1), ] <- c("Power Law distribution", info)

info <- c(round(coef(cities_spl, "sigma")[2], 3), round(AIC(cities_spl), 3))
results[(nrow(results) +  1), ] <- c("Shifted Power Law distribution", info)

results
