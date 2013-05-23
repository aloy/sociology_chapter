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
autism.full$sicdegp <- relevel(autism.full$sicdegp, ref = 3)

# Center age at 2 years to interpret the intercept
autism.full$age2 <- autism.full$age - 2

# Taking the complete cases for our example
autism.modframe <- subset(autism.full, 
	select = c(childid, sicdegp, age2, vsae, gender, race, bestest2))
autism.modframe <- na.omit(autism.modframe)

#-------------------------------------------------------------------------------
# Model selection/checking: Lineups vs. conventional inference
#-------------------------------------------------------------------------------

### Initial model
(M1 <- lmer(vsae ~ age2 + (age2 - 1 | childid), data = autism.modframe))


### Is the linear random component for age enough?

## Conventional answer
M2 <- lmer(vsae ~ age2 + (age2 + I(age2^2) - 1 | childid), data = autism.modframe)
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

### Do we need to allow for correlation between the two random effects?

## Conventional tests say yes, but less convincingly
M3 <- lmer(vsae ~ age2 + (age2 - 1 | childid) + (I(age2^2) - 1 | childid), data = autism.modframe)
anova(M2, M3)

