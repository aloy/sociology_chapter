#-------------------------------------------------------------------------------
# Script generating figure 7
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

library(checkpoint)
checkpoint("2016-06-10")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(grid)

# Calculating the no. of observations for each county
radon <- ddply(radon, .(county), transform, n = length(county))

# The model using only counties with at least 4 obs.
fm <- lmer(log.radon ~ basement + uranium + (basement | county), 
           data = subset(radon, n > 4))

#-------------------------------------------------------------------------------
# Cyclone plots
#-------------------------------------------------------------------------------

# The *observed residuals* and other variables needed for lineup
observed.cyclone.df <- data.frame(fm@frame, resid = resid(fm))

# Simulating a set of 19 null plots
sim.y   <- simulate(fm, nsim = 19, seed = 987654321) ## Change seed for diff. nulls                       
sim.mod <- apply(sim.y, 2, refit, object = fm)

# Extract level-1 residuals both from the simulated models
null.cyclone.df <- ldply(sim.mod, function(x) data.frame(x@frame, resid = resid(x)))
null.cyclone.df$.id <- as.numeric( str_extract(null.cyclone.df$.id, "\\d+") )

# Combine for plotting
cyclone.df <- rbind.fill(observed.cyclone.df, null.cyclone.df)

# Randomizing the panels
plot.order <- sample.int(20, 20)
location <- plot.order[1] # location of the true plot
cyclone.df$.id <- rep(plot.order, each = nrow(fm@frame))

# Ordering the subjects by IQR
cyclone.df <- ddply(cyclone.df, .(.id, county), transform, iqr = IQR(resid, na.rm=TRUE))
cyclone.df <- ddply(cyclone.df, .(.id), transform, rank = order(order(iqr, county)))
cyclone.df <- ddply(cyclone.df, .(.id, county), transform, rank = min(rank))
cyclone.df <- dlply(cyclone.df, .(.id), transform, rank = rank(rank))
cyclone.df <- ldply(cyclone.df, function(df) {
  df$rank <- factor(df$rank)
  levels(df$rank) <- 1:length(levels(df$rank))
  return(df)
})

# Rendering the lineup
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