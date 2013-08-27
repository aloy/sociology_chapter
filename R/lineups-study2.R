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

load("dialyzer-heterogeneous2.RData")

lineup.df <- rbind.fill(m2.sim.resids, m2.resid.df)
plot.order <- sample.int(20, 20)
lineup.df$.n <- rep(plot.order, each = 140)

location <- plot.order[20]
qplot(x = QB, y = resid, data = lineup.df, geom = "boxplot", facets = ~ .n, 
      fill = QB, outlier.size = 1.5, alpha=I(0.6)) + 
	xlab(NULL) + 
	ylab(NULL) + 
	scale_fill_brewer("", palette="Set2") +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


##### lineup number 4 #####

load("dialyzer-heterogeneous2.RData")

lineup.df <- rbind.fill(m2.sim.resids, m2.resid.df)
plot.order <- sample.int(20, 20)
lineup.df$.n <- rep(plot.order, each = 140)

location <- plot.order[20]
qplot(x = QB, y = resid, data = lineup.df, geom = "jitter", facets = ~ .n, 
      colour = QB, outlier.size = 1.5, alpha=I(0.6)) + 
	xlab(NULL) + 
	ylab(NULL) + 
	scale_colour_brewer("", palette="Set2") +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


##### lineup number 5 #####

# Function to construct the envelopes typically used with Q-Q plots
env <- function(x, conf = .95){
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

# This is a reprise of Q-Q plots from the study Lendie and Heike ran,
# but this time we are framing it in the same way as the other lineups
# for this paper---i.e., we have five sets of nulls.

# This loads all 5 sets of nulls
load("radon-normal-qq-sim.RData")

lineup.df <- sim_slope_qq_df_1 # simply change the last number here (1-5)
lineup.df <- subset(lineup.df, select = -std.slope)
plot.order <- sample.int(20, 20)
lineup.df$.n <- rep(plot.order, each = 85)

location <- plot.order[20]

# adding the envelope info
lineup.df <- ddply(lineup.df[complete.cases(lineup.df),], .(.n), 
                   transform, band = env(`slope`), 
                   x = qqnorm(`slope`, plot.it=FALSE)$x,
                   x.env = sort(qqnorm(`slope`, plot.it=FALSE)$x))

# adding the reference line
lineup.df <- ddply(lineup.df, .(.n), transform,
                   y.line = HLMdiag:::qqlineInfo(slope)[1] + HLMdiag:::qqlineInfo(slope)[2] * x.env)


ggplot(lineup.df, aes(x = x, y = slope)) + 
	facet_wrap(~ .n, ncol = 5) + 
	geom_ribbon(aes(x = lineup.df$x.env, ymin = lineup.df $band.lower, 
	                ymax = lineup.df$band.upper), alpha = .2) + 
	geom_line(aes(x = x.env, y = y.line), colour = I("grey60")) + 
	geom_point() + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
	      axis.title.x = element_blank(), axis.title.y = element_blank())


##### lineup number 6 #####

# I think it would be useful to add in standardized versions of the
# Q-Q plots from lineup 5. This would help add some support for what
# Heike and I were working on directly before my defense.

# This loads all 5 sets of nulls
load("radon-normal-qq-sim.RData")

lineup.df <- sim_slope_qq_df_1 # simply change the last number here (1-5)
plot.order <- sample.int(20, 20)
lineup.df$.n <- rep(plot.order, each = 85)

location <- plot.order[20]

# adding the envelope info
lineup.df <- ddply(lineup.df[complete.cases(lineup.df),], .(.n), 
                   transform, band = env(std.slope), 
                   x = qqnorm(std.slope , plot.it=FALSE)$x,
                   x.env = sort(qqnorm(std.slope, plot.it=FALSE)$x))

# adding the reference line
lineup.df <- ddply(lineup.df, .(.n), transform,
                   y.line = HLMdiag:::qqlineInfo(std.slope)[1] + 
                   HLMdiag:::qqlineInfo(std.slope)[2] * x.env)

ggplot(lineup.df, aes(x = x, y = std.slope)) + 
	facet_wrap(~ .n, ncol = 5) + 
	geom_ribbon(aes(x = lineup.df$x.env, ymin = lineup.df $band.lower, 
	                ymax = lineup.df$band.upper), alpha = .2) + 
	geom_line(aes(x = x.env, y = y.line), colour = I("grey60")) + 
	geom_point() + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
	      axis.title.x = element_blank(), axis.title.y = element_blank())


##### lineup number 7 #####

# If we look at standardized random effects for the normal setting,
# then we should rerun the situation where the random effects
# folow a t(3) distribution using standardized random effects.