#-------------------------------------------------------------------------------
# Script for constructing the lineups for the second MTurk study for
# "Visual Inference for Linear Mixed-Effects Models"
#
# Adam Loy, Heike Hofmann, Dianne Cook
# June 2013
#-------------------------------------------------------------------------------

library(grid)
library(ggplot2)
library(plyr)

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

location <- plot.order[1]
qplot(x = factor(rank), y = EB.resid, data = resids, geom = "boxplot", outlier.size = 1.5) + 
	coord_flip() + 
	ylim(-, ) + 
	xlab(NULL) + 
	ylab(NULL) + 
	theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
	      axis.text.x = element_blank(), axis.ticks.y = element_blank(), 
	      panel.grid.major.y = element_blank()) +
	facet_wrap(~ sample)

