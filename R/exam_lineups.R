#-------------------------------------------------------------------------------
# Script constructing lineups for use in a paper demonstrating their utility
# in checking the assumptions for hierarchical linear models.
#
# Adam Loy
# May 2013
#
# EXAM DATA
# Data description:
# 4059 students nested within 65 schools
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

setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")
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

#-------------------------------------------------------------------------------
# Lineup to test for random slope
#-------------------------------------------------------------------------------
# Initial model
M1 <- lmer(normexam ~ standLRT + (1 | school), data = Exam)

# Do we need a random slope? - We can compare the observed data to simulated data
# to see if the model can generate data that look like the original.
qplot(x = standLRT, y = normexam, data = Exam, geom = "smooth", group = school, se = F, method = "lm", colour = I("black"))

# The "easy way" with smoothers, but alpha refers to the bands...
m1.sims  <- simulate(M1, nsim = 19, seed = 1234)
m1.refit <- lapply(m1.sims, refit, object = M1)
m1.simy <- lapply(m1.refit, function(x) x@y)

sim.y <- do.call("cbind", m1.simy)
sim.y <- melt(sim.y)[,-1]
names(sim.y) <- c(".n", "y")
sim.y$.n <- as.numeric(str_extract(sim.y$.n, "\\d+"))
sim.y$standLRT <- rep(Exam$standLRT, rep = 19)
sim.y$school <- rep(Exam$school, rep = 19)

true.y <- data.frame(y = Exam$normexam, standLRT = Exam$standLRT, school = Exam$school)

# qplot(x = standLRT, y = y, data = true.y, group = school, geom = "smooth", 
#	method = "lm", se=F, colour = I("black"), alpha = I(0.5)) %+% 
#	lineup(true = true.y, samples = sim.y) + 
#	facet_wrap( ~ .sample, ncol=5) + 
#	xlab("standardized LRT score") + 
#	ylab("GCSE exam score")
	
# The "harder way" to get alpha working
m1.fitted <- ddply(sim.y, .(.n, school), function(x) {
	m <- lm(y ~ standLRT, data = x)
	data.frame(x, fitted = fitted(m))
})

true.df <- subset(Exam, select = c(school, normexam, standLRT))
colnames(true.df)[2] <- "y"
true.fitted <- ddply(true.df, .(school), function(x) {
	m <- lm(y ~ standLRT, data = x)
	data.frame(x, fitted = fitted(m))	
})

qplot(x = standLRT, y = fitted, data = true.fitted, group = school, 
	geom = "line", alpha = I(0.5)) %+% 
	lineup(true = true.fitted, samples = m1.fitted) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
		axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#	xlab("standardized LRT score") + 
#	ylab("GCSE exam score")


## Adding the random slope
M2 <- lmer(normexam ~ standLRT + (standLRT  | school), data = Exam)

M2.sims  <- simulate(M2, nsim = 19, seed = 1234)
M2.refit <- lapply(M2.sims, refit, object = M2)
M2.simy <- lapply(M2.refit, function(x) x@y)

sim2.y <- do.call("cbind", M2.simy)
sim2.y <- melt(sim2.y)[,-1]
names(sim2.y) <- c(".n", "y")
sim2.y$.n <- as.numeric(str_extract(sim2.y$.n, "\\d+"))
sim2.y$standLRT <- rep(Exam$standLRT, rep = 19)
sim2.y$school <- rep(Exam$school, rep = 19)

# The "harder way" to get alpha working
M2.fitted <- ddply(sim2.y, .(.n, school), function(x) {
	m <- lm(y ~ standLRT, data = x)
	data.frame(x, fitted = fitted(m))
})

true.df <- subset(Exam, select = c(school, normexam, standLRT))
colnames(true.df)[2] <- "y"

true.fitted <- ddply(true.df, .(school), function(x) {
	m <- lm(y ~ standLRT, data = x)
	data.frame(x, fitted = fitted(m))	
})

# save(M2.fitted, true.fitted, file = "exam-fanned-withslope.RData")

qplot(x = standLRT, y = fitted, data = true.fitted, group = school, 
	geom = "line", alpha = I(0.5)) %+% 
	lineup(true = true.fitted, samples = M2.fitted) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
		axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


#-------------------------------------------------------------------------------
# Lineups to test for need for correlated random effects
#-------------------------------------------------------------------------------

# Do we need to allow for correlation? Typically you would add the random effect, then
# test for the correlation. So that is what we will try.
M3 <- lmer(normexam ~ standLRT + (standLRT - 1 | school) + (1 | school), data = Exam)

## The lineup -need to compare simulated ranefs to ranefs of M2
set.seed(987654321)
M3.sims  <- simulate(M3, nsim = 19)
M3.refit <- lapply(M3.sims, refit, object = M3)
M3.sim.ranef <- lapply(M3.refit, function(x) ranef(x)[[1]])

M3.sim.ranef <- do.call("rbind", M3.sim.ranef)
M3.sim.ranef$.n <- rownames(M3.sim.ranef)
M3.sim.ranef$.n <- as.numeric(str_extract(M3.sim.ranef$.n, "\\d+"))

true.M2.ranef <- ranef(M2)$school

qplot(x = `(Intercept)`, y = standLRT, data = true.M2.ranef, 
	geom = c("point", "smooth"), method = "lm", se = F, alpha = I(0.4)) %+% 
	lineup(true = true.M2.ranef, samples = M3.sim.ranef) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
ggsave("exam_corr_lineup7.pdf")
save(true.M2.ranef, M3.sim.ranef, file="exam-corr.RData")

location <- 7
make_interactive(filename= sprintf("exam-corr-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("exam-corr-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")


## A simulation when we know that we DO NOT need correlation
set.seed(9221632)
M3.sims2  <- simulate(M3, nsim = 20)
M3.refit2 <- lapply(M3.sims2[,-20], refit, object = M3)
M3.sim.ranef2 <- lapply(M3.refit2, function(x) ranef(x)[[1]])

M3.sim.ranef2 <- do.call("rbind", M3.sim.ranef2)
M3.sim.ranef2$.n <- rownames(M3.sim.ranef2)
M3.sim.ranef2$.n <- as.numeric(str_extract(M3.sim.ranef2$.n, "\\d+"))

uncorr.mod <- refit(M2, M3.sims2[,20])
true.uncorr.ranef <- ranef(uncorr.mod)$school

qplot(x = `(Intercept)`, y = standLRT, data = true.uncorr.ranef, 
	geom = c("point", "smooth"), method = "lm", se = F, alpha = I(0.4)) %+% 
	lineup(true = true.uncorr.ranef, samples = M3.sim.ranef2) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

#-------------------------------------------------------------------------------
# Lineup to test for level-1 heteroscedasticity
#-------------------------------------------------------------------------------
### Negative
fm2 <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + 
              (standLRT | school), Exam, REML = FALSE)

sim_fm2 <- simulate(fm2, nsim = 19)
refit_fm2 <- apply(sim_fm2, 2, refit, object = fm2)

## Extract level-1 residuals
sim_fm2_lev1_resid <- ldply(refit_fm2, function(x) 
	data.frame(x@frame, resid = resid(x)))

## Labeling for nullabor
sim_fm2_lev1_resid$.n <- rep(1:19, each = 4059)
names(sim_fm2_lev1_resid)[4:5] <- c("standLRT2", "standLRT3")

## Formatting for nullabor
lev1_resid_fm2 <- data.frame(fm2@frame, resid = resid(fm2))
names(lev1_resid_fm2)[3:4] <- c("standLRT2", "standLRT3")
class(lev1_resid_fm2[,3])  <- "numeric"

## Creating lineup
qplot(standLRT2, resid, data = lev1_resid_fm2,
      	geom = "point", alpha = I(0.4)) %+%
  	lineup(true = lev1_resid_fm2, samples = sim_fm2_lev1_resid) +
  	facet_wrap(~ .sample, ncol = 5) +
  	geom_hline(aes(yintercept = 0), colour = I("white"), alpha = I(0.5)) + 
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

location <- 5
make_interactive(filename= sprintf("exam-homogeneity-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("exam-homogeneity-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")
