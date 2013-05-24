#-------------------------------------------------------------------------------
# Script constructing lineups for use in a paper demonstrating their utility
# in checking the assumptions for hierarchical linear models.
#
# Adam Loy
# May 2013
#
# EXAM DATA
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

setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")
library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)

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

qplot(x = standLRT, y = y, data = true.y, group = school, geom = "smooth", 
	method = "lm", se=F, colour = I("black"), alpha = I(0.5)) %+% 
	lineup(true = true.y, samples = sim.y) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab("standardized LRT score") + 
	ylab("GCSE exam score")
	
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
	xlab("standardized LRT score") + 
	ylab("GCSE exam score")


# Do we need to allow for correlation? Typically you would add the random effect, then
# test for the correlation. So that is what we will try.
M2 <- lmer(normexam ~ standLRT + (standLRT  | school), data = Exam)
M3 <- lmer(normexam ~ standLRT + (standLRT - 1 | school) + (1 | school), data = Exam)

## The lineup -need to compare simulated ranefs to ranefs of M2
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
	xlab("age - 2") + 
	ylab(expression(paste((age - 2)^2)))
