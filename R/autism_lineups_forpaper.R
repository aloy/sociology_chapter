-------------------------------------------------------------------------------
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

# Center age at 2 years to interpret the intercept
autism.full$age2 <- autism.full$age - 2

# Taking the complete cases for our example
autism.modframe <- subset(autism.full, 
	select = c(childid, sicdegp, age2, vsae, gender, race, bestest2))
autism.modframe <- na.omit(autism.modframe)

#-------------------------------------------------------------------------------
# Model selection Lineups vs. conventional inference
#-------------------------------------------------------------------------------

set.seed(9221632)

M1 <- lmer(vsae ~ age2 + ( age2 - 1 | childid ), data = autism.modframe)

### sicdegp
# We can look at the level-1 residuals, though they may be less interesting
M1.sim.sicdegp  <- data.lineup.explvar1(null.model = M1, variable = "sicdegp", 
	data = autism.modframe, std = FALSE)
M1.sim.sicdegp.std  <- data.lineup.explvar1(null.model = M1, variable = "sicdegp",
	 data = autism.modframe, std = TRUE)
M1.true.sicdegp <- data.frame(residual = resid(M1), 
	sicdegp = autism.modframe$sicdegp)
M1.true.sicdegp.std <- data.frame(residual = HLMresid(M1, level = 1, 
	standardize = TRUE), sicdegp = autism.modframe$sicdegp)


qplot(x = sicdegp, y = residual, data = M1.true.sicdegp, geom = "boxplot", 
	fill = sicdegp, outlier.size = 0) %+% 
	lineup(true = M1.true.sicdegp, samples = M1.sim.sicdegp) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10) + 
	ylab(NULL) + xlab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
	legend.position = "none")


# mixing up the order
autism.modframe$sicdegp <- relevel(autism.modframe$sicdegp, ref = 3)

### sicdegp
# We can look at the level-1 residuals, though they may be less interesting
M1.sim.sicdegp  <- data.lineup.explvar1(null.model = M1, variable = "sicdegp", 
	data = autism.modframe, std = FALSE)
M1.sim.sicdegp.std  <- data.lineup.explvar1(null.model = M1, variable = "sicdegp",
	 data = autism.modframe, std = TRUE)
M1.true.sicdegp <- data.frame(residual = resid(M1), 
	sicdegp = autism.modframe$sicdegp)
M1.true.sicdegp.std <- data.frame(residual = HLMresid(M1, level = 1, 
	standardize = TRUE), sicdegp = autism.modframe$sicdegp)


qplot(x = sicdegp, y = residual, data = M1.true.sicdegp, geom = "boxplot", 
	fill = sicdegp, outlier.size = 0) %+% 
	lineup(true = M1.true.sicdegp, samples = M1.sim.sicdegp) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10) + 
	ylab(NULL) + xlab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
	legend.position = "none")
