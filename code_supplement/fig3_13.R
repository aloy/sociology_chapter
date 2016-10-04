#-------------------------------------------------------------------------------
# Script generating figures 3 and 12
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------
library(checkpoint)
checkpoint("2016-06-10")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals and the data set
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)

# Function used to simulate from the null model and extract the level-1 residuals
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

#-------------------------------------------------------------------------------
# Figure 3
#-------------------------------------------------------------------------------

# The null model
M2 <- lmer(vsae ~ poly(age2, 2) + (age2 -1 | childid) + (I(age2^2) - 1 | childid), 
           data = autism)

# Simulating the level-1 residuals for the null plots
set.seed(9221632)
M2.sim.sicdegp  <- data.lineup.explvar1(null.model = M2, variable = "sicdegp", 
                                        data = autism, std = FALSE)

# Extracting the level-1 residuals for the true plot
M2.true.sicdegp <- data.frame(residual = resid(M2), sicdegp = autism$sicdegp)

qplot(x = sicdegp, y = residual, data = M2.true.sicdegp, geom = "boxplot", 
      fill = sicdegp, outlier.size = 2, alpha=I(0.6)) %+% 
  lineup(true = M2.true.sicdegp, samples = M2.sim.sicdegp) + 
  facet_wrap( ~ .sample, ncol=5) + 
  ylim(-10, 10) + 
  ylab(NULL) + 
  xlab(NULL) + 
  scale_fill_brewer("", palette="Set2", labels=c("low", "medium", "high")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


#-------------------------------------------------------------------------------
# Figure 12
#-------------------------------------------------------------------------------

# Model with only a linear random slope
M1 <- lmer(vsae ~ poly(age2, 2) + (age2 - 1 | childid), data = autism)

# Model with a quadratic random slope
M2 <- lmer(vsae ~ poly(age2, 2) + (age2 -1 | childid) + (I(age2^2) - 1 | childid), 
           data = autism)


# Simulating data for the null plots from M1
M1.sims  <- simulate(M1, nsim = 19, seed = 12345)
M1.refit <- lapply(M1.sims, refit, object = M1)
M1.sim.y <- lapply(M1.refit, function(x) getME(x, "y"))

# Formatting the simulated data for use in a lineup
M1.sim.y <- do.call("cbind", M1.sim.y)
M1.sim.y <- melt(M1.sim.y)[,-1]
names(M1.sim.y) <- c(".n", "y")
M1.sim.y$y[M1.sim.y$y < 0] <- 0
# M1.sim.y$vsae <- rep(autism$vsae, rep = 19)
M1.sim.y$childid <- rep(autism$childid, rep = 19)
M1.sim.y$age2 <- rep(autism$age2, rep = 19)
M1.sim.y$.n <- as.numeric(str_extract(M1.sim.y$.n, "\\d+"))

# Formatting the true data
autism.true.y <- data.frame(y = autism$vsae, age2 = autism$age2, childid = autism$childid)

# Adding the true data to the data frame for the null plots
true_pos <- sample(20, 1) # print to reveal position
fig12_df <- nullabor:::add_true(M1.sim.y, autism.true.y, pos = true_pos)

# Creating the lineup
ggplot(data = fig12_df, aes(x = age2, y = y, group = childid)) +
  geom_line(alpha = 0.4) +
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())