#-------------------------------------------------------------------------------
# Script constructing lineups of the Q-Q plots.
#
# Adam Loy
# May 2013
#
# RADON DATA
# 
#
# References: - Gelman, A & Pardoe, I (2006). Bayesian measures of explained
#				variance and pooling in multilevel (hierarchical) models.
#				Technometrics, 48(2), 241--251
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")
library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(LearnBayes)

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

sim_t_hlm <- function(.mod) {
	vc <- VarCorr( .mod )
	D  <- as.matrix( bdiag(vc) )
	sig.e <- sigma(.mod)
	
	dims <- .mod@dims
	n <- dims[["n"]]
	m <- dims[["q"]] / nrow(D)

	## normal errors
	e  <- rnorm(n = n, mean = 0, sd = sig.e)

	## mutlivariate t random effects
	b <- rmt(n = m, mean = c(0, 0), S = D, df = 3)
	
	## Generating y
	bvec <- c(b[,1], b[,2])
	y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% bvec + e
	
	return( as.numeric(y) )
}


radon <- read.csv("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter/data/radon_for_sims.csv")

fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

#-------------------------------------------------------------------------------
# Q-Q plots -- observed radon data
#-------------------------------------------------------------------------------
b <- ranef(fm)[[1]] # notice that this is actually a matrix

sim.y   <- simulate(fm, nsim = 19)                        
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

sim.b0 <- llply(sim.mod, function(x) ranef(x)[[1]][,1])   ## a list of random slopes
sim.b0 <- melt( do.call("rbind", sim.b0) )[,-2]           ## changing to a data frame
names(sim.b0) <- c("sample", "(Intercept)")                 
sim.b0        <- arrange(sim.b0, sample)                  ## ordering by simulation

sim.b0$.n <- as.numeric( str_extract(sim.b0$sample, "\\d+") )
sim.b0 <- ddply(sim.b0, .(.n), transform, band = sim_env(`(Intercept)`), 
	x = sort(qqnorm(`(Intercept)`, plot.it=FALSE)$x))

b0 <- transform(b, band = sim_env(`(Intercept)`), 
	x = sort(qqnorm(`(Intercept)`, plot.it=FALSE)$x))

# Lineup of random intercepts
qplot(sample = X.Intercept., data = b0, stat = "qq") %+%
	lineup(true = b0, sample = sim.b0) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles")



sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

b1 <- transform(b, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

# Lineup of random slopes
qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles")


#-------------------------------------------------------------------------------
# Q-Q plots -- simulating normal random effects
#-------------------------------------------------------------------------------

b <- ranef(fm)[[1]] # notice that this is actually a matrix

sim.y   <- simulate(fm, nsim = 20)                        
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

sim.b0 <- llply(sim.mod, function(x) ranef(x)[[1]][,1])   ## a list of random slopes
sim.b0 <- melt( do.call("rbind", sim.b0) )[,-2]           ## changing to a data frame
names(sim.b0) <- c("sample", "(Intercept)")                 
sim.b0        <- arrange(sim.b0, sample)                  ## ordering by simulation

sim.b0$.n <- as.numeric( str_extract(sim.b0$sample, "\\d+") )
sim.b0 <- ddply(sim.b0, .(.n), transform, band = sim_env(`(Intercept)`), 
	x = sort(qqnorm(`(Intercept)`, plot.it=FALSE)$x))

b0 <- subset(sim.b0, .n == 20)
sim.b0 <- sim.b0[-which(sim.b0$.n == 20),] 

# Lineup of random intercepts
qplot(sample = X.Intercept., data = b0, stat = "qq") %+%
	lineup(true = b0, sample = sim.b0) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles")



sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

b1 <- subset(sim.b1, .n == 20)
sim.b1 <- sim.b1[-which(sim.b1$.n == 20),] 

# Lineup of random slopes
qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles")


#-------------------------------------------------------------------------------
# Q-Q plots -- simulating non-normal random effects
#-------------------------------------------------------------------------------

y.b.t <- sim_t_hlm(fm)
refit.b.t <-  refit(fm, y.b.t)
b.t <- ranef(refit.b.t)[[1]]

sim.y   <- simulate(fm, nsim = 19)                        
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

sim.b0 <- llply(sim.mod, function(x) ranef(x)[[1]][,1])   ## a list of random slopes
sim.b0 <- melt( do.call("rbind", sim.b0) )[,-2]           ## changing to a data frame
names(sim.b0) <- c("sample", "(Intercept)")                 
sim.b0        <- arrange(sim.b0, sample)                  ## ordering by simulation

sim.b0$.n <- as.numeric( str_extract(sim.b0$sample, "\\d+") )
sim.b0 <- ddply(sim.b0, .(.n), transform, band = sim_env(`(Intercept)`), 
	x = sort(qqnorm(`(Intercept)`, plot.it=FALSE)$x))

b0 <- transform(b.t, band = sim_env(`(Intercept)`), 
	x = sort(qqnorm(`(Intercept)`, plot.it=FALSE)$x))


qplot(sample = X.Intercept., data = b0, stat = "qq") %+%
	lineup(true = b0, sample = sim.b0) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme(panel.margin = unit(0, "lines"))



sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

b1 <- transform(b.t, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

# Lineup of random slopes
qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles")
