# Analysis of the patent data set, please download the file here:
# www.stat.uni-muenchen.de/~kneib/regressionsbuch/download/patentdata.raw

# Prepare data

patents <- read.table("patentdata.raw", header = TRUE, sep = "\t")

patents$biopharm <- as.factor(patents$biopharm)
patents$ustwin <- as.factor(patents$ustwin)
patents$patus <- as.factor(patents$patus)
patents$patgsgr <- as.factor(patents$patgsgr)

# Plot histograms on a linear scale and on a log-log scale. The red line in
# the right plot seems a pretty good fit, thus 10 could be the lower bound
# of the Power Law behavior

breaks <- (max(patents$nclaims)-min(patents$nclaims))
hist <- hist(patents$nclaims, breaks, plot = FALSE)

par(mfrow = c(1, 2))

plot(
  x    = hist$mids[1:100],
  y    = hist$counts[1:100],
  type = "l",
  main = "Linear Scale",
  xlab = "Patent Claims",
  ylab = "Count",
)

plot(
  x    = hist$mids[1:100],
  y    = hist$counts[1:100],
  type = "l",
  main = "Log-Log Scale",
  xlab = "Patent Claims",
  ylab = "Count",
  log  = "xy"
)

x <- hist$mids[10:100][hist$counts[10:100] != 0]
y <- hist$counts[10:100][hist$counts[10:100] != 0]
coef <- lm(log10(y) ~ log10(x))$coefficients

abline(coef, col = "red", lwd = 2)

# Look for optimal lower bound minimizing the KS-statistic. According to the
# KS-statistic, 31 is the optimal lower bound. The bootstrap test indictates
# that we cannot reject the Power Law hypothesis

library(poweRlaw)

m <- displ(patents$nclaims)
est <- estimate_xmin(m)

m$setXmin(10)
set.seed(1337)
bootstrap_p(m)

m$setXmin(est)
set.seed(1337)
bootstrap_p(m)

# Perform two standard regression with the data, first with nclaims as response
# variable and then with log(nclaims) as response. The first model violates the
# normal distribution assumption for the error term. However, the logarithmic
# model yields good results

patents_normal <- lm(
  formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  data    = patents
)

plot(density(patents$nclaims, bw = 1), main = "nclaims")
plot(resid(patents_normal), main = "nclaims as response", ylab = "Residuals")

AIC(patents_normal)
BIC(patents_normal)

patents_log <- lm(
  formula = log(nclaims) ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  data    = patents
)

AIC(patents_log)
BIC(patents_log)

# Examine the data set with GAMLSS and compare different skew distributions.
# The log-normal distribution performs best, but the Shifted Power Law also
# works out well

library(gamlss)

source("../distributions/PL.R")
source("../distributions/SPL.R")
source("../misc/gamlss.summary.R")

patents_exp <- gamlss(
  formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family  = EXP,
  data    = patents,
  control = gamlss.control(n.cyc = Inf)
)

AIC(patents_exp)
BIC(patents_exp)

patents_lno <- gamlss(
  formula       = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = LNO,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_lno)
BIC(patents_lno)

patents_pareto2o <- gamlss(
  formula       = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = PARETO2o,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_pareto2o)
BIC(patents_pareto2o)

patents_pl <- gamlss(
  formula       = nclaims ~ 1,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = PL,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_pl)
BIC(patents_pl)

patents_spl <- gamlss(
  formula       = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = SPL,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_spl)
BIC(patents_spl)

patents_wei <- gamlss(
  formula       = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = WEI,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_wei)
BIC(patents_wei)

# Drop the explanatory variables for the lower bound of the Shifted Power Law
# model. This way, we avoid problems with the interpretation

patents_spl_alpha <- gamlss(
  formula       = nclaims ~ 1,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = SPL,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_spl_alpha)
BIC(patents_spl_alpha)

# Estimate two GAMLSS with a lower bound of 31. Without a lower bound, GAMLSS
# with the non-shifted Power Law distribution do not perform well

patents <- patents[patents$nclaims >= 31, ]

patents_pl_bounded <- gamlss(
  formula       = nclaims ~ 1,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = PL,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_pl_bounded)
BIC(patents_pl_bounded)

patents_spl_bounded <- gamlss(
  formula       = nclaims ~ 1,
  sigma.formula = nclaims ~ biopharm + ustwin + patus + patgsgr + year + ncountry,
  family        = SPL,
  data          = patents,
  control       = gamlss.control(n.cyc = Inf)
)

AIC(patents_spl_bounded)
BIC(patents_spl_bounded)
