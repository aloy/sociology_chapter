#-------------------------------------------------------------------------------
# Script generating figures 6, 9, 12, 14, 15
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------
library(checkpoint)
checkpoint("2016-06-10")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals and data set
library(MEMSS)    # for the Dialyzer data set
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)


#-------------------------------------------------------------------------------
# Figure 6: Lineup to test for homogeneity
#-------------------------------------------------------------------------------

# Fitting the quartic LME model
M2 <- lmer(rate ~ (pressure + I(pressure^2) + I(pressure^3) + I(pressure^4))*QB + 
             (pressure + I(pressure^2) | Subject), data = Dialyzer)

# Extracting the residuals
m2.resid.df <- data.frame(M2@frame, resid = resid(M2))

# Simulating residuals for the null plots from model M2
m2.sims <- simulate(M2, nsim = 19)
m2.refit <- lapply(m2.sims, refit, object = M2)
m2.sim.resids <- ldply(m2.refit, function(x) data.frame(x@frame, resid = resid(x)))
m2.sim.resids$.n <- as.numeric(str_extract(m2.sim.resids$.id, "\\d+"))

# Creating the lineup
qplot(pressure, resid, data = m2.resid.df,
      geom = "point", alpha = I(0.5)) %+%
  lineup(true = m2.resid.df, samples = m2.sim.resids) +
  facet_wrap(~ .sample, ncol = 5) +
  xlab(NULL) + 
  ylab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

#-------------------------------------------------------------------------------
# Figure 9: Lineup to test for linearity
#-------------------------------------------------------------------------------

# Fitting the quadratic model
M1 <- lmer(rate ~ (pressure + I(pressure^2))*QB + (pressure + I(pressure^2) | Subject), 
           data = Dialyzer)

# Extracting the residuals
m1.resid.df <- data.frame(M1@frame, resid = resid(M1))

# Simulating residuals for the null plots from model M1
set.seed("20040110")
m1.sims <- simulate(M1, nsim = 19)
m1.refit <- lapply(m1.sims, refit, object = M1)
m1.sim.resids <- ldply(m1.refit, function(x) data.frame(x@frame, resid = resid(x)))
m1.sim.resids$.n <- as.numeric(str_extract(m1.sim.resids$.id, "\\d+"))

resids <- rbind.fill(m1.resid.df, m1.sim.resids)

# Randomizing the panels
plot.order <- sample.int(20, 20)
location <- plot.order[1] # location of the true plot
resids$.id <- rep(plot.order, each = nrow(M1@frame))


# Creating the lineup
ggplot(resids, aes(x = pressure, y = resid)) + 
  geom_point() + 
  geom_smooth(method = "loess") +
  facet_wrap(~ .id, ncol = 5) +
  xlab(NULL) + 
  ylab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

#-------------------------------------------------------------------------------
# Figure 12: Lineup to test for homogeneity
#-------------------------------------------------------------------------------

# Discretizing pressure for boxplots representation
fig11_df <- m2.sim.resids
fig11_df$pressurecat <- cut(fig11_df$pressure, 
                         breaks = c(0, 0.45, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25))

fig11_true <- m2.resid.df
fig11_true$pressurecat <- cut(fig11_true$pressure, breaks = c(0, 0.45, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25))

# Adding the true data to the data frame for the null plots
true_pos <- sample(20, 1) # print to reveal position
fig11_df <- nullabor:::add_true(fig11_df, fig11_true, pos = true_pos)

# Creating the lineup
qplot(x = pressurecat, y = resid, data = fig11_df, geom = "boxplot", fill = pressurecat, alpha=I(0.8)) +
  facet_wrap(~ .sample, ncol = 5) +
  xlab(NULL) + 
  ylab(NULL) + 
  scale_fill_brewer() +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.position="none") + 
  scale_x_discrete(expand = c(.1,.1))


#-------------------------------------------------------------------------------
# Figure 14: Lineup to test for homogeneity
#-------------------------------------------------------------------------------

# Formatting data frame for plotting and sorting Subjects by residual variance
true_pos <- sample(20, 1) # print to reveal position
fig13_df <- nullabor:::add_true(m1.sim.resids, m1.resid.df, pos = true_pos)

fig13_df <- ddply(fig13_df, .(.sample, Subject), transform, var = var(resid))
fig13_df <- ddply(fig13_df, .(.sample), transform, subid = rank(var, ties.method="min"))


qplot(subid, resid, data = fig13_df, geom = c("point"), fill=I("grey70"), alpha=I(0.5)) +
  facet_wrap(~ .sample, ncol = 5) +
  geom_hline(yintercept=0, colour="red") +
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.position="none") +
  scale_x_continuous(expand = c(.075,.075))

#-------------------------------------------------------------------------------
# Figure 15: Lineup to test for homogeneity
#-------------------------------------------------------------------------------

qplot(Subject, resid, data = m1.resid.df, geom = "boxplot", fill=I("grey70"), alpha=I(0.8)) %+%
  lineup(true = m1.resid.df, samples = m1.sim.resids) +
  facet_wrap(~ .sample, ncol = 5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position="none") +
  scale_x_discrete(expand = c(.075,.075))