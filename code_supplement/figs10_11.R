#-------------------------------------------------------------------------------
# Script generating figures 10 and 11
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

library(checkpoint)
checkpoint("2016-06-10")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(mvtnorm)

# Function simulating the confidence envelope used for Q-Q plots
sim_env <- function(x, conf = .95){
  n <- length(x)
  P <- ppoints(x)
  z <- qnorm(P)
  a <- as.numeric(HLMdiag:::qqlineInfo(x)[1])
  b <- as.numeric(HLMdiag:::qqlineInfo(x)[2])
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- a + b * z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
  return(data.frame(lower, upper))
}

# Function simulating data with mutlivariate t random effects and normal errors
sim_t_hlm <- function(.mod) {
  vc <- VarCorr( .mod )
  D  <- as.matrix( bdiag(vc) )
  sig.e <- sigma(.mod)
  
  n <- getME(.mod, "n")
  m <- getME(.mod, "q") / nrow(D)
  
  ## normal errors
  e  <- rnorm(n = n, mean = 0, sd = sig.e)
  
  ## mutlivariate t random effects
  b <- rmvt(n = m, sigma = D, df = 3)
  
  ## Generating y
  bvec <- c(b[,1], b[,2])
  y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% bvec + e
  
  return( as.numeric(y) )
}


# Load the data
data(radon, "HLMdiag")

#-------------------------------------------------------------------------------
# Figure 10: Lineup of Q-Q plots for the random slope
#-------------------------------------------------------------------------------

# Fit the proposed model
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

# Extract the random slope
b <- ranef(fm)[[1]] # notice that this is actually a matrix

# Simulating data from the model
sim.y   <- simulate(fm, nsim = 19, seed = 987654321)                        
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

# Extracting the random slopes from the simulations
# and formating a data set of the 20 sets of Q-Q plots
sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1[complete.cases(sim.b1),], .(.n), transform, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

b1 <- transform(b, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

# Rendering the lineup
qplot(sample = basement, data = b1, stat = "qq") %+%
  lineup(true = b1, sample = sim.b1) + 
  facet_wrap(~ .sample, ncol = 5) + 
  geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
  ylab(NULL) + 
  xlab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

#-------------------------------------------------------------------------------
# Figure 11: Lineup of Q-Q plots for the random slope simulated from a t dsn
#-------------------------------------------------------------------------------

# Fitting the proposed model
fm <- lmer(log.radon ~ basement + uranium + (1 | county) + 
             (basement - 1 | county), data = radon)

# Simulating data with random effects from a t distribution
set.seed(-2029298609)
y.b.t <- sim_t_hlm(fm)
refit.b.t <-  refit(fm, y.b.t)
b.t <- ranef(refit.b.t)[[1]]

# Simulating a normal LME model
sim.y   <- simulate(fm, nsim = 19)                        
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

# Extracting the random slopes from the simulations
# and formating a data set of the 20 sets of Q-Q plots
sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

b1 <- transform(b.t, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

# Rendering the lineup of the random slope
location <- 6
qplot(sample = basement, data = b1, stat = "qq") %+%
  lineup(true = b1, sample = sim.b1, pos=location) + 
  facet_wrap(~ .sample, ncol = 5) + 
  geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
  #	xlab("Normal Quantiles") + ylab("Sample Quantiles") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())