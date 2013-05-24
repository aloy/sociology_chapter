#-------------------------------------------------------------------------------
# Script constructing lineups for use in a paper demonstrating their utility
# in checking the assumptions for hierarchical linear models.
#
# Adam Loy
# May 2013
#
# Data description: 
# prospective longitudinal study
# 158 children who were classified as autism spectrum disorder children.
# Variables: 
# * childid - child id
# * vsae - Vineland socialization age equivalent
# * sicdegp - Sequenced inventory of communication development expressive group,
#	i.e., expressive language group (1 = low, 2 = medium, 3 = high)
# * age - ages (2, 3, 5, 9, 13)
# * agemom - age of mother when child was born
# * gender - gender (1 = male, 2 = female)
# * race - race (1 = white, 2 = non-white)
# * bestest2 - diagnosis at age 2 (1 = autism, 2 = PDD, 3 = non-spectrum)
# * bestest9 - diagnosis at age 9 (1 = autism, 2 = PDD, 3 = non-spectrum)
# * bestviq - verbal IQ at age 2
# * bestnviq - non-verbal IQ at age 2
#
# Sources:
# West, Welch, and Galecki (2006).
#   Linear Mixed Models: A Practical Guide Using Statistical Software.
#
# Anderson, D. K., Oti, R. S., Lord, C., & Welch, K. (2009).
#	Patterns of Growth in Adaptive Social Abilities Among Children with 
#	Autism Spectrum Disorders. Journal of Abnormal Child Psychology, 37(7), 
#	1019–1034. doi:10.1007/s10802-009-9326-0
# 
# Anderson, D. K., Lord, C., Risi, S., DiLavore, P. S., Shulman, C., 
#	Thurm, A., et al. (2007). Patterns of growth in verbal abilities among 
#	children with autism spectrum disorder. Journal of Consulting and Clinical 
#	Psychology, 75(4), 594–604. doi:10.1037/0022-006X.75.4.594
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")
library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(WWGbook)  # for half of data set - see ?autism
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)

# Function to make lineups for explanatory variables less tedious - level 1
data.lineup.explvar1 <- function(null.model, variable, data, nsim = 19, std = FALSE) {
	mod.sims  <- simulate(null.model, nsim = nsim)
	mod.refit <- lapply(mod.sims, refit, object = null.model)
	if(std){
		mod.sim.resid <- lapply(mod.refit, HLMresid, level = 1, 
			type = "EB", standardize = TRUE)
	} else{
		mod.sim.resid <- lapply(mod.refit, resid)
	}

	mod.sim.resid <- do.call("cbind", mod.sim.resid)
	mod.sim.resid <- melt(mod.sim.resid)[,-1]
	names(mod.sim.resid) <- c(".n", "residual")
	mod.sim.resid$.n <- as.numeric(str_extract(mod.sim.resid $.n, "\\d+"))
	
	mod.sim.resid <- cbind(mod.sim.resid, data[, variable])
	names(mod.sim.resid)[-c(1:2)] <- variable
	
	return(mod.sim.resid)
}

# Function to make lineups for explanatory variables less tedious - level 2
data.lineup.explvar2 <- function(null.model, variable, data, nsim = 19) {
	mod.sims  <- simulate(null.model, nsim = nsim)
	mod.refit <- lapply(mod.sims, refit, object = null.model)
	mod.sim.resid <- lapply(mod.refit, HLMresid, level = "childid")

	mod.sim.resid <- do.call("rbind", mod.sim.resid)
	mod.sim.resid$.n <- as.numeric(str_extract(rownames(mod.sim.resid), "\\d+"))
	mod.sim.resid$childid <- str_extract(rownames(mod.sim.resid), "[.]\\d+")
	mod.sim.resid$childid <- as.numeric(str_extract(mod.sim.resid$childid, "\\d+"))
	
	mod.sim.resid <- cbind(mod.sim.resid, data[, variable])
	names(mod.sim.resid)[-c(1:4)] <- variable
	
	return(mod.sim.resid)
}


# Reading in more of the autism data set
demographics <- read.csv("data/autism_demographics.csv")

# combining the data sets
autism.full <- merge(x = autism, y = demographics, by.x = c("childid", "sicdegp"),
	by.y = c("newinit", "sicdegp"), all.x = TRUE)
	
# Creating factors where necessary
autism.full$sicdegp  <- factor(autism.full$sicdegp, labels = c("low", "med", "high"))
autism.full$gender   <- factor(autism.full$gender, labels = c("male", "female"))
autism.full$race     <- factor(autism.full$race, labels = c("white", "nonwhite"))
autism.full$bestest2 <- factor(autism.full$bestest2, labels = c("autism", "pdd"))

# Revel sicdegp
# autism.full$sicdegp <- relevel(autism.full$sicdegp, ref = 3)

# Center age at 2 years to interpret the intercept
autism.full$age2 <- autism.full$age - 2

# Taking the complete cases for our example
autism.modframe <- subset(autism.full, 
	select = c(childid, sicdegp, age2, vsae, gender, race, bestest2))
autism.modframe <- na.omit(autism.modframe)

#-------------------------------------------------------------------------------
# Model selection Lineups vs. conventional inference
#-------------------------------------------------------------------------------

### Selecting the random effects structure
#-------------------------------------------------------------------------------

### Initial model
(M1 <- lmer(vsae ~ poly(age2, 2) + (age2 - 1 | childid), data = autism.modframe))


### Is the linear random component for age enough?

## Conventional answer
M2 <- lmer(vsae ~ poly(age2, 2) + (age2 -1 | childid) + (I(age2^2) - 1 | childid), 
	data = autism.modframe)
anova(M1, M2)

## Lineup comparing the observed growth curves to versions simulated under M1
M1.sims  <- simulate(M1, nsim = 19, seed = 12345)
M1.refit <- lapply(M1.sims, refit, object = M1)
M1.sim.y <- lapply(M1.refit, function(x) x@y)

M1.sim.y <- do.call("cbind", M1.sim.y)
M1.sim.y <- melt(M1.sim.y)[,-1]
names(M1.sim.y) <- c(".n", "y")
M1.sim.y$y[M1.sim.y$y < 0] <- 0
M1.sim.y$.n <- as.numeric(str_extract(M1.sim.y$.n, "\\d+"))
M1.sim.y$vsae <- rep(autism.modframe$vsae, rep = 19)
M1.sim.y$childid <- rep(autism.modframe$childid, rep = 19)
M1.sim.y$age2 <- rep(autism.modframe$age2, rep = 19)

autism.true.y <- data.frame(y = autism.modframe$vsae, age2 = autism.modframe$age2, 
	childid = autism.modframe$childid)

qplot(x = age2, y = y, data = autism.true.y, group = childid, 
		geom = "line", se=F, alpha = I(0.3)) %+% 
	lineup(true = autism.true.y, samples = M1.sim.y) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylab("VSAE") + 
	xlab("age - 2")
	
### Now we can check the same type of lineup to see if this seems sufficient
M2.sims  <- simulate(M2, nsim = 19)
M2.refit <- lapply(M2.sims, refit, object = M2)
M2.sim.y <- lapply(M2.refit, function(x) x@y)

M2.sim.y <- do.call("cbind", M2.sim.y)
M2.sim.y <- melt(M2.sim.y)[,-1]
names(M2.sim.y) <- c(".n", "y")
M2.sim.y$y[M2.sim.y$y < 0] <- 0
M2.sim.y$.n <- as.numeric(str_extract(M2.sim.y$.n, "\\d+"))
M2.sim.y$vsae <- rep(autism.modframe$vsae, rep = 19)
M2.sim.y$childid <- rep(autism.modframe$childid, rep = 19)
M2.sim.y$age2 <- rep(autism.modframe$age2, rep = 19)

qplot(x = age2, y = y, data = autism.true.y, group = childid, 
		geom = "line", se=F, alpha = I(0.3)) %+% 
	lineup(true = autism.true.y, samples = M2.sim.y) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylab("VSAE") + 
	xlab("age - 2")


### Do we need to allow for correlation between the two random effects?

## Conventional tests say yes, but less convincingly
M3 <- lmer(vsae ~ age2 + I(age2^2) + (age2 - 1 | childid) + (I(age2^2) - 1 | childid), data = autism.modframe)

# M3 <- lmer(vsae ~ age2 + I(age2^2) + (age2 + I(age2^2) - 1 | childid), data = autism.modframe)

anova(M2, M3)

## The lineup -need to compare simulated ranefs to ranefs of M2
M3.sims  <- simulate(M3, nsim = 19)
M3.refit <- lapply(M3.sims, refit, object = M3)
M3.sim.ranef <- lapply(M3.refit, function(x) ranef(x)[[1]])

M3.sim.ranef <- do.call("rbind", M3.sim.ranef)
M3.sim.ranef$.n <- rownames(M3.sim.ranef)
M3.sim.ranef$.n <- as.numeric(str_extract(M3.sim.ranef$.n, "\\d+"))

true.M2.ranef <- ranef(M2)$childid 

qplot(x = age2, y = `I(age2^2)`, data = true.M2.ranef, 
	geom = c("point", "smooth"), method = "lm", se = F, alpha = I(0.4)) %+% 
	lineup(true = true.M2.ranef, samples = M3.sim.ranef) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab("age - 2") + 
	ylab(expression(paste((age - 2)^2)))


### Selecting the fixed effects structure
#-------------------------------------------------------------------------------

### Let's assume that we had not already hypothesize a quadratic relationship
### between VSAE and age. Then we may have started with the following model:

M1.alt <- lmer(vsae ~ age2 + ( age2 - 1 | childid ), data = autism.modframe)

M1.alt.sim.age  <- data.lineup.explvar1(null.model = M1.alt, variable = "age2", 
	data = autism.modframe, std = FALSE)
M1.alt.true.sim.age <- data.frame(residual = HLMresid(M2, level = 1), 
	age2 = autism.modframe$age2)


qplot(x = age2, y = residual, data = M1.alt.true.sim.age, 
	geom = c("jitter", "smooth")) %+% 
	lineup(true = M1.alt.true.sim.age, samples = M1.alt.sim.age) + 
	facet_wrap( ~ .sample, ncol=5)


### Level-2 variables (child-level)

child.df <- subset(autism.modframe, 
	select = c(childid, sicdegp, gender, race, bestest2))
child.df <- unique(child.df)

### sicdegp
M2.ml <- update(M2, . ~ ., REML = FALSE)
M4 <- lmer(vsae ~ age2 + I(age2^2) + sicdegp + (age2 + I(age2^2) - 1 | childid), data = autism.modframe, REML = FALSE)
anova(M2.ml, M4)

# We can look at the level-1 residuals, though they may be less interesting
M2.sim.sicdegp  <- data.lineup.explvar1(null.model = M2, variable = "sicdegp", 
	data = autism.modframe, std = FALSE)
M2.sim.sicdegp.std  <- data.lineup.explvar1(null.model = M2, variable = "sicdegp",
	 data = autism.modframe, std = TRUE)
M2.true.sicdegp <- data.frame(residual = resid(M2), 
	sicdegp = autism.modframe$sicdegp)
M2.true.sicdegp.std <- data.frame(residual = HLMresid(M2, level = 1, 
	standardize = TRUE), sicdegp = autism.modframe$sicdegp)

qplot(x = sicdegp, y = residual, data = M2.true.sicdegp, geom = "boxplot", 
	fill = sicdegp, outlier.size = 0) %+% 
	lineup(true = M2.true.sicdegp, samples = M2.sim.sicdegp) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10) + 
	ylab("level-1 residuals")

qplot(x = sicdegp, y = residual, data = M2.true.sicdegp.std, geom = "boxplot", 
	fill = sicdegp, outlier.size = 0) %+% 
	lineup(true = M2.true.sicdegp.std, samples = M2.sim.sicdegp.std) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-4, 4) + 
	ylab("standardized residual")
	
# The level-2 residuals
M2.sim.sicdegp2  <- data.lineup.explvar2(null.model = M2, variable = "sicdegp", 
	data = child.df)
M2.true.sicdegp2 <- data.frame(ranef(M2)[[1]], sicdegp = child.df$sicdegp)
names(M2.true.sicdegp2)[2] <- "I(age2^2)"

qplot(x = sicdegp, y = age2, data = M2.true.sicdegp2, geom = "boxplot", 
	fill = sicdegp, outlier.size = 0) %+% 
	lineup(true = M2.true.sicdegp2, samples = M2.sim.sicdegp2) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10)

qplot(x = sicdegp, y = `I(age2^2)`, data = M2.true.sicdegp2, geom = "boxplot", 
	fill = sicdegp, outlier.size = 0) %+% 
	lineup(true = M2.true.sicdegp2, samples = M2.sim.sicdegp2) + 
	facet_wrap( ~ .sample, ncol=5)


### gender
M5 <- lmer(vsae ~ age2 + I(age2^2) + sicdegp + gender + (age2 + I(age2^2) - 1 | childid), data = autism.modframe, REML = FALSE)
anova(M4, M5)

M4.sim.gender  <- data.lineup.explvar1(null.model = M4, variable = "gender", 
	data = autism.modframe, std = FALSE)
M4.sim.gender.std  <- data.lineup.explvar1(null.model = M4, variable = "gender",
	 data = autism.modframe, std = TRUE)
M4.true.gender <- data.frame(residual = resid(M4), 
	gender = autism.modframe$gender)
M4.true.gender.std <- data.frame(residual = HLMresid(M4, level = 1, 
	standardize = TRUE), gender = autism.modframe$gender)

qplot(x = gender, y = residual, data = M4.true.gender, geom = "boxplot", 
	fill = gender, outlier.size = 0) %+% 
	lineup(true = M4.true.gender, samples = M4.sim.gender) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10)

qplot(x = gender, y = residual, data = M4.true.gender.std, geom = "boxplot", 
	fill = gender, outlier.size = 0) %+% 
	lineup(true = M4.true.gender.std, samples = M4.sim.gender.std) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-4, 4) + 
	ylab("standardized residual")

# The level-2 residuals
M2.sim.gender2  <- data.lineup.explvar2(null.model = M2, variable = "gender", 
	data = child.df)
M2.true.gender2 <- data.frame(ranef(M2)[[1]], gender = child.df$gender)
names(M2.true.gender2)[2] <- "I(age2^2)"

qplot(x = gender, y = age2, data = M2.true.gender2, geom = "boxplot", 
	fill = gender, outlier.size = 0) %+% 
	lineup(true = M2.true.gender2, samples = M2.sim.gender2) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10)

qplot(x = gender, y = `I(age2^2)`, data = M2.true.gender2, geom = "boxplot", 
	fill = gender, outlier.size = 0) %+% 
	lineup(true = M2.true.gender2, samples = M2.sim.gender2) + 
	facet_wrap( ~ .sample, ncol=5)


### race
M6 <- lmer(vsae ~ age2 + I(age2^2) + sicdegp + race + (age2 + I(age2^2) - 1 | childid), data = autism.modframe, REML = FALSE)
anova(M4, M6)

M4.sim.race  <- data.lineup.explvar1(null.model = M4, variable = "race", 
	data = autism.modframe, std = FALSE)
M4.sim.race.std  <- data.lineup.explvar1(null.model = M4, variable = "race",
	 data = autism.modframe, std = TRUE)
M4.true.race <- data.frame(residual = resid(M4), 
	race = autism.modframe$race)
M4.true.race.std <- data.frame(residual = HLMresid(M4, level = 1, 
	standardize = TRUE), race = autism.modframe$race)


qplot(x = race, y = residual, data = M4.true.race.std, geom = "boxplot", 
	fill = race, outlier.size = 0) %+% 
	lineup(true = M4.true.race.std, samples = M4.sim.race.std) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-5, 5) + 
	ylab("standardized residual")

# The level-2 residuals
M2.sim.race2  <- data.lineup.explvar2(null.model = M2, variable = "race", 
	data = child.df)
M2.true.race2 <- data.frame(ranef(M2)[[1]], race = child.df$race)
names(M2.true.race2)[2] <- "I(age2^2)"

qplot(x = race, y = age2, data = M2.true.race2, geom = "boxplot", 
	fill = race, outlier.size = 0) %+% 
	lineup(true = M2.true.race2, samples = M2.sim.race2) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10)

qplot(x = race, y = `I(age2^2)`, data = M2.true.race2, geom = "boxplot", 
	fill = race, outlier.size = 0) %+% 
	lineup(true = M2.true.race2, samples = M2.sim.race2) + 
	facet_wrap( ~ .sample, ncol=5)

### bestest2
M7 <- lmer(vsae ~ age2 + I(age2^2) + sicdegp + bestest2 + (age2 + I(age2^2) - 1 | childid), data = autism.modframe, REML = FALSE)
anova(M4, M7)

M4.sim.bestest2  <- data.lineup.explvar1(null.model = M4, variable = "bestest2", 
	data = autism.modframe, std = FALSE)
M4.sim.bestest2.std  <- data.lineup.explvar1(null.model = M4, variable = "bestest2",
	 data = autism.modframe, std = TRUE)
M4.true.bestest2 <- data.frame(residual = resid(M4), 
	bestest2 = autism.modframe$bestest2)
M4.true.bestest2.std <- data.frame(residual = HLMresid(M4, level = 1, 
	standardize = TRUE), bestest2 = autism.modframe$bestest2)

qplot(x = bestest2, y = residual, data = M4.true.bestest2.std, geom = "boxplot", 
	fill = bestest2, outlier.size = 0) %+% 
	lineup(true = M4.true.bestest2.std, samples = M4.sim.bestest2.std) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-5, 5) + 
	ylab("standardized residual")
	
# The level-2 residuals
M2.sim.bestest22  <- data.lineup.explvar2(null.model = M2, variable = "bestest2", 
	data = child.df)
M2.true.bestest22 <- data.frame(ranef(M2)[[1]], bestest2 = child.df$bestest2)
names(M2.true.bestest22)[2] <- "I(age2^2)"

qplot(x = bestest2, y = age2, data = M2.true.bestest22, geom = "boxplot", 
	fill = bestest2, outlier.size = 0) %+% 
	lineup(true = M2.true.bestest22, samples = M2.sim.bestest22) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10)

qplot(x = bestest2, y = `I(age2^2)`, data = M2.true.bestest22, geom = "boxplot", 
	fill = bestest2, outlier.size = 0) %+% 
	lineup(true = M2.true.bestest22, samples = M2.sim.bestest22) + 
	facet_wrap( ~ .sample, ncol=5)

#-------------------------------------------------------------------------------
# Model checking
#-------------------------------------------------------------------------------

### Normality of the random effects

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

radon <- read.csv("~/Documents/Thesis/Dissertation/eresids-chapter/minconfounded_chapter/data/radon_for_sims.csv") ## find radon_for_sims.csv

fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

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

b1 <- transform(b, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme(panel.margin = unit(0, "lines"))


D <- bdiag( VarCorr(fm) )














M2.sim.ranef <- lapply(M2.refit, function(x) ranef(x)[[1]]) # M2.refit from above

M2.sim.ranef <- do.call("rbind", M2.sim.ranef)
M2.sim.ranef$.n <- as.numeric(str_extract(rownames(M2.sim.ranef), "\\d+"))
M2.sim.ranef$childid <- str_extract(rownames(M2.sim.ranef), "[.]\\d+")
M2.sim.ranef$childid <- as.numeric(str_extract(M2.sim.ranef$childid, "\\d+"))

M2.true.ranef <- ranef(M2)[[1]]

qplot(sample = age2, data = M2.true.ranef, stat = "qq") %+%
	lineup(true = M2.true.ranef, sample = M2.sim.ranef) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
	xlab("Normal Quantiles") + ylab("Sample Quantiles") +  
	theme_bw() + 
	theme(panel.margin = unit(0, "lines"))


### Within-group homoscedasticity
qplot(x = factor(childid), y = resid(M4), geom = "boxplot", data = autism.modframe)
