#-------------------------------------------------------------------------------
# Script generating figure 7
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------
library(checkpoint)
checkpoint("2016-06-10")

library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(ggplot2)
library(grid)
library(plyr)
library(MASS)

# Fitting the proposed model
model <- lmer(sbvalue ~ treatment + week*baseline + I(week^2)*baseline + 
                (0 + week + I(week^2) | subject), data = ahd)

#-------------------------------------------------------------------------------
# Cyclone plots
#-------------------------------------------------------------------------------

# Checking between group homogeneity
## Residuals from true model
set.seed(123456789)

## Simulating null data and refitting the model
nsim <- 19
refit_fm <- vector(mode = "list", length = nsim)
for(i in 1:nsim) {
  refit_fm[[i]] <- refit(model, simulate(model))
}

## Extract level-1 residuals both from the sims and the true model
sim_resids <- ldply(refit_fm, function(x){
  e <- resid(x)
  subj <- x@frame$subject
  return(data.frame(EB.resid = e, subject = subj))
})
sim_resids <- data.frame(.id = rep(1:19, each = nrow(model@frame)), sim_resids)

true_resids <- data.frame(.id = 20, EB.resid = resid(model), subject = model@frame$subject)

# combining the results for plotting
resids <- rbind(sim_resids, true_resids)

# randomizing the order of the plots within the lineup
resids$sample <- sample(20, 20, replace = FALSE)[resids$.id]

# Ordering the subjects by IQR
resids <- ddply(resids, .(.id, subject), transform, iqr = IQR(EB.resid))
resids <- ddply(resids, .(.id), transform, rank = order(order(iqr, subject)))
resids <- ddply(resids, .(.id, subject), transform, rank = min(rank))
resids <- ddply(resids, .(.id), transform, rank = rank(rank))

# Rendering the lineup
qplot(x = factor(rank), y = EB.resid, data = resids, 
      geom = "boxplot", xlab = "", ylab = "", outlier.size = 1.5) + 
  coord_flip() + 
  ylim(-125, 125) + 
  theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.y = element_blank()) +
  facet_wrap(~sample)