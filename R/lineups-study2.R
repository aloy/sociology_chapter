#-------------------------------------------------------------------------------
# Script for constructing more lineups for the second MTurk study for
# "Visual Inference for Linear Mixed-Effects Models"
#
# Adam Loy, Heike Hofmann, Dianne Cook
# August 2013
#
# This script gives code constructing lineups for the second Mturk study
# that will help "fill out" the paper.
#-------------------------------------------------------------------------------

library(grid)
library(ggplot2)
library(plyr)
library(nullabor)

##### lineup number 1 #####

load("exam-fanned-withslope.RData")

location <- 18
qplot(x = standLRT, y = fitted, data = true.fitted, group = school, 
	geom = "line", alpha = I(0.5)) %+% 
	lineup(true = true.fitted, samples = M2.fitted, pos = 18) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
		axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


##### lineup number 2 #####

load("radon-cyclone.RData")

cyclone.df <- rbind.fill(observed.cyclone.df, null.cyclone.df)

plot.order <- sample.int(20, 20)
cyclone.df$.id <- rep(plot.order, each = 919)

cyclone.df <- ddply(cyclone.df, .(.id, county), transform, iqr = IQR(resid))
cyclone.df <- ddply(cyclone.df, .(.id), transform, rank = order(order(iqr, county)))
cyclone.df <- ddply(cyclone.df, .(.id, county), transform, rank = min(rank))
cyclone.df <- dlply(cyclone.df, .(.id), transform, rank = rank(rank))
cyclone.df <- ldply(cyclone.df, function(df) {
	df$rank <- factor(df$rank)
	levels(df$rank) <- 1:length(levels(df$rank))
	return(df)
})

location <- plot.order[1]
qplot(x = factor(rank), y = resid, data = cyclone.df, geom = "boxplot", outlier.size = 1.5) + 
	coord_flip() + 
	ylim(-2, 2) + 
	xlab(NULL) + 
	ylab(NULL) + 
	facet_wrap(~ .id) + 
	theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
	      axis.text.x = element_blank(), axis.ticks.y = element_blank(), 
	      panel.grid.major.y = element_blank())


##### lineup number 3 #####

