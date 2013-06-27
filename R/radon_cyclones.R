#-------------------------------------------------------------------------------
# Script for constructing the lineups of cyclone plots for the radon data
#
# Adam Loy, Heike Hofmann, Dianne Cook
# June 2013
#-------------------------------------------------------------------------------

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(grid)

radon <- ddply(radon, .(county), transform, n = length(county))

# The model - with fewer counties
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = subset(radon, n > 4))

# The *observed residuals* and other variables needed for lineup
observed.cyclone.df <- data.frame(fm@frame, resid = resid(fm))

# Simulating a set of 19 null plots
sim.y   <- simulate(fm, nsim = 19, seed = 987654321) ## Change seed for diff. nulls                       
sim.mod <- apply(sim.y, 2, refit, object = fm)

null.cyclone.df <- ldply(sim.mod, function(x) data.frame(x@frame, resid = resid(x)))
null.cyclone.df$.id <- as.numeric( str_extract(null.cyclone.df$.id, "\\d+") )

cyclone.df <- rbind.fill(observed.cyclone.df, null.cyclone.df)

plot.order <- sample.int(20, 20)
location <- plot.order[1]
cyclone.df$.id <- rep(plot.order, each = nrow(fm@frame))

cyclone.df <- ddply(cyclone.df, .(.id, county), transform, iqr = IQR(resid))
cyclone.df <- ddply(cyclone.df, .(.id), transform, rank = order(order(iqr, county)))
cyclone.df <- ddply(cyclone.df, .(.id, county), transform, rank = min(rank))
cyclone.df <- dlply(cyclone.df, .(.id), transform, rank = rank(rank))
cyclone.df <- ldply(cyclone.df, function(df) {
	df$rank <- factor(df$rank)
	levels(df$rank) <- 1:length(levels(df$rank))
	return(df)
})

qplot(x = factor(rank), y = resid, data = cyclone.df, geom = "boxplot", outlier.size = 1.5) + 
	coord_flip() + 
	ylim(-2, 2) + 
	xlab(NULL) + 
	ylab(NULL) + 
	theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
	      axis.text.x = element_blank(), axis.ticks.y = element_blank(),
	      axis.ticks.x = element_blank(), 
	      panel.grid.major.y = element_blank()) +
	facet_wrap(~ .id)


save(observed.cyclone.df, null.cyclone.df, file = "radon-cyclone.RData")