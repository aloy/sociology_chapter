#-------------------------------------------------------------------------------
# Script generating figures 1, 3, and 4
#
# Testing for necessity of a random slope
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------
library(checkpoint)
checkpoint("2016-06-10")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)

#-------------------------------------------------------------------------------
# Figure 1
#-------------------------------------------------------------------------------
## Initial model
M1 <- lmer(normexam ~ standLRT + (1 | school), data = Exam)

# Do we need a random slope? - We can compare the observed data to simulated data
# to see if the model can generate data that look like the original.
qplot(x = standLRT, y = normexam, data = Exam, geom = "smooth", group = school, 
      se = F, method = "lm", colour = I("black"))

# Simulating data for the null plots from model M1
m1.sims  <- simulate(M1, nsim = 19, seed = 1234)
m1.refit <- lapply(m1.sims, refit, object = M1)
m1.simy <- lapply(m1.refit, function(x) getME(x, "y"))

# Formatting a data frame for nullabor
sim.y <- do.call("cbind", m1.simy)
sim.y <- melt(sim.y)[,-1]
names(sim.y) <- c(".n", "y")
sim.y$.n <- as.numeric(str_extract(sim.y$.n, "\\d+"))
sim.y$standLRT <- rep(Exam$standLRT, rep = 19)
sim.y$school <- rep(Exam$school, rep = 19)

true.y <- data.frame(y = Exam$normexam, standLRT = Exam$standLRT, school = Exam$school)

# Extracting information for the linear smoothers than will be plotted for each group
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

# Creating the lineup: linear smoothers are plotted for each group,
# the different panes represent different data sets.
qplot(x = standLRT, y = fitted, data = true.fitted, group = school, 
      geom = "line", alpha = I(0.5)) %+% 
  lineup(true = true.fitted, samples = m1.fitted) + 
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#	xlab("standardized LRT score") + 
#	ylab("GCSE exam score")

#-------------------------------------------------------------------------------
# Figure 3
#-------------------------------------------------------------------------------

## Adding the random slope and refitting the LME model
M2 <- lmer(normexam ~ standLRT + (standLRT  | school), data = Exam)

# Simulating data for the null plots from model M1
M2.sims  <- simulate(M2, nsim = 19, seed = 1234)
M2.refit <- lapply(M2.sims, refit, object = M2)
M2.simy <- lapply(M2.refit, function(x) getME(M2, "y"))

sim2.y <- do.call("cbind", M2.simy)
sim2.y <- melt(sim2.y)[,-1]
names(sim2.y) <- c(".n", "y")
sim2.y$.n <- as.numeric(str_extract(sim2.y$.n, "\\d+"))
sim2.y$standLRT <- rep(Exam$standLRT, rep = 19)
sim2.y$school <- rep(Exam$school, rep = 19)

# Extracting information for the linear smoothers than will be plotted for each group
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

qplot(x = standLRT, y = fitted, data = true.fitted, group = school, 
      geom = "line", alpha = I(0.5)) %+% 
  lineup(true = true.fitted, samples = M2.fitted) + 
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

#-------------------------------------------------------------------------------
# Figure 4: Testing the need for correlated random effects
#-------------------------------------------------------------------------------
# Fitting a model with independent random effects
M3 <- lmer(normexam ~ standLRT + (standLRT - 1 | school) + (1 | school), data = Exam)

# Simulating random effects from model M3 with indepenent random effects
set.seed(987654321)
M3.sims  <- simulate(M3, nsim = 19)
M3.refit <- lapply(M3.sims, refit, object = M3)
M3.sim.ranef <- lapply(M3.refit, function(x) ranef(x)[[1]])

# Formatting a data frame of the null plots for plotting
M3.sim.ranef <- do.call("rbind", M3.sim.ranef)
M3.sim.ranef$.n <- rownames(M3.sim.ranef)
M3.sim.ranef$.n <- as.numeric(str_extract(M3.sim.ranef$.n, "\\d+"))

# Extracting the random effects from the model where we allow random 
# effects to be correlated
true.M2.ranef <- ranef(M2)$school

# Adding the true data to the data frame for the null plots
true_pos <- sample(20, 1) # print to reveal position
fig4_df <- nullabor:::add_true(M3.sim.ranef, true.M2.ranef, pos = true_pos)

# Creating the lineup
ggplot(data = fig4_df, aes(x = `(Intercept)`, y = standLRT)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.4) +
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())