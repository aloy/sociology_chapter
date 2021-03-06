#-------------------------------------------------------------------------------
# Script constructing lineups for use in a paper demonstrating their utility
# in checking the assumptions for hierarchical linear models.
#
# Adam Loy
# May 2013
#
# DIALYZER DATA
# Description: 4059 students nested within 65 schools
# Variables - 
# *school - school id
# *student - student id
# *normexam - student's standardized exam score at 16
# *schgend - school's gender
# *schavg - 
# *vr -
# *intake -
# *standLRT - student's standardized score on the London reading test (age 11)
# *sex - student's gender
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

setwd("~/Documents/Research/sociology_chapter/")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)

source("R/add_interaction.R")

library(MEMSS) # for the Dialyzer data set

#-------------------------------------------------------------------------------
# Lineup to test for linearity
#-------------------------------------------------------------------------------

M1 <- lmer(rate ~ (pressure + I(pressure^2))*QB + (pressure + I(pressure^2) | Subject), 
	data = Dialyzer)

m1.resid.df <- data.frame(pressure = M1@frame$pressure, resid = resid(M1))

m1.sims <- simulate(M1, nsim = 19)
m1.refit <- lapply(m1.sims, refit, object = M1)
m1.sim.resids <- ldply(m1.refit, function(x) data.frame(x@frame, resid = resid(x)))
m1.sim.resids$.n <- as.numeric(str_extract(m1.sim.resids$.id, "\\d+"))

location <- sample(20, 1)
qplot(pressure, resid, data = m1.resid.df, geom = c("point", "smooth")) %+%
  	lineup(true = m1.resid.df, samples = m1.sim.resids, pos=location) +
  	facet_wrap(~ .sample, ncol = 5) +
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
ggsave(file=sprintf("dialyzer-nonlinear-%s.pdf", location))
make_interactive(filename= sprintf("dialyzer-nonlinear-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("dialyzer-nonlinear-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")

# saving data
# 
save(m1.resid.df, m1.sim.resids, file = "dialyzer-nonlinear.RData")

#-------------------------------------------------------------------------------
# Lineup to test for homogeneity
#-------------------------------------------------------------------------------

# Model on page 217 of Pinheiro and Bates

M2 <- lmer(rate ~ (pressure + I(pressure^2) + I(pressure^3) + I(pressure^4))*QB + (pressure + I(pressure^2) | Subject), data = Dialyzer)

m2.resid.df <- data.frame(M2@frame, resid = resid(M2))

m2.sims <- simulate(M2, nsim = 19)
m2.refit <- lapply(m2.sims, refit, object = M2)
m2.sim.resids <- ldply(m2.refit, function(x) data.frame(x@frame, resid = resid(x)))
m2.sim.resids$.n <- as.numeric(str_extract(m2.sim.resids$.id, "\\d+"))

qplot(pressure, resid, data = m2.resid.df,
      	geom = "point", alpha = I(0.5)) %+%
  	lineup(true = m2.resid.df, samples = m2.sim.resids) +
  	facet_wrap(~ .sample, ncol = 5) +
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

location <- 3
make_interactive(filename= sprintf("dialyzer-heterogeneous-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("dialyzer-heterogeneous-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")


# saving data
# save(m2.resid.df, m2.sim.resids, file = file.choose(new = TRUE))


### Another look at homogeneity using a different explanatory variable


# Lineup using boxplots
lineup.df <- rbind.fill(m2.sim.resids, m2.resid.df)
plot.order <- sample.int(20, 20)
lineup.df$.n <- rep(plot.order, each = 140)

qplot(x = QB, y = resid, data = lineup.df, geom = "boxplot", facets = ~ .n, fill = QB,
		outlier.size = 1.5, alpha=I(0.6)) + 
	xlab(NULL) + 
	ylab(NULL) + 
	scale_fill_brewer("", palette="Set2") +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

# Lineup using dotplots
lineup.df <- rbind.fill(m2.sim.resids, m2.resid.df)
plot.order <- sample.int(20, 20)
lineup.df$.n <- rep(plot.order, each = 140)

qplot(x = QB, y = resid, data = lineup.df, geom = "jitter", facets = ~ .n, 
      colour = QB, outlier.size = 1.5, alpha=I(0.6)) + 
	xlab(NULL) + 
	ylab(NULL) + 
	scale_colour_brewer("", palette="Set2") +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
