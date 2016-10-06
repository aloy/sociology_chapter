#-------------------------------------------------------------------------------
# Script generating figures and tables in the article
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------
library(checkpoint) # For reproducibility
checkpoint("2016-06-10")

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)
library(mvtnorm)
library(parallel)
require(magrittr)
require(xtable)
# devtools::install_github("heike/vinference")
library(vinference) 

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


# Function simulating the confidence envelope used for Q-Q plots
sim_env <- function(x, conf = .95){
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

# Function simulating data with mutlivariate t random effects and normal errors
sim_t_hlm <- function(.mod) {
  vc <- VarCorr( .mod )
  D  <- as.matrix( bdiag(vc) )
  sig.e <- sigma(.mod)
  
  n <- getME(.mod, "n")
  m <- getME(.mod, "q") / nrow(D)
  
  ## normal errors
  e  <- rnorm(n = n, mean = 0, sd = sig.e)
  
  ## mutlivariate t random effects
  b <- rmvt(n = m, sigma = D, df = 3)
  
  ## Generating y
  bvec <- c(b[,1], b[,2])
  y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% bvec + e
  
  return( as.numeric(y) )
}

# Function calculating the residual std. dev. for separate LM fits to each group
calc_s_df <- function(x, formula) {
  mod <- lm(formula, data = x)
  smry <- summary(mod)
  s  <- smry$sigma  # residual std. dev.
  df <- smry$df[2]  # n_i - r_i
  
  return(c(s = s, df = df))
}

# Function calculating the test statistic, df, and p-value for standard test
het.chisq.test <- function(cutoff) {
  # Model from which cyclone lineup was created
  fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = subset(radon, n >= cutoff))
  
  # Calculating the d_i from equation (6)
  test.df2 <- ddply(fm@frame, .(county), calc_s_df, formula = log.radon ~ basement)
  test.df2 <- transform(test.df2, d = ( log(s^2) - ( sum(df * log(s^2)) / sum(df) ) ) / sqrt( 2 / df ) )
  
  # Calculating the test stat from equation (7) along with the standard p-value
  H2 <- sum(test.df2$d^2)
  df <- nrow(test.df2) - 1
  pval <- pchisq(H2, df = df, lower.tail = FALSE)
  
  return(c(H = H2, df = df, p.val = pval))
}



#-------------------------------------------------------------------------------
# Lineups based on the Exam data set
# Figures 2, 4 and 5
#-------------------------------------------------------------------------------

## Initial model
exam.mod1 <- lmer(normexam ~ standLRT + (1 | school), data = Exam)

# Do we need a random slope? - We can compare the observed data to simulated data
# to see if the model can generate data that look like the original.
qplot(x = standLRT, y = normexam, data = Exam, geom = "smooth", group = school, 
      se = F, method = "lm", colour = I("black"))

# Simulating data for the null plots from model M1
exam.m1.sims  <- simulate(exam.mod1, nsim = 19, seed = 1234)
exam.m1.refit <- lapply(exam.m1.sims, refit, object = exam.mod1)
exam.m1.simy <- lapply(exam.m1.refit, function(x) getME(x, "y"))

# Formatting a data frame for nullabor
exam.sim.y <- do.call("cbind", exam.m1.simy)
exam.sim.y <- melt(exam.sim.y)[,-1]
names(exam.sim.y) <- c(".n", "y")
exam.sim.y$.n <- as.numeric(str_extract(exam.sim.y$.n, "\\d+"))
exam.sim.y$standLRT <- rep(Exam$standLRT, rep = 19)
exam.sim.y$school <- rep(Exam$school, rep = 19)

exam.true.y <- data.frame(y = Exam$normexam, standLRT = Exam$standLRT, school = Exam$school)

# Extracting information for the linear smoothers than will be plotted for each group
exam.m1.fitted <- ddply(exam.sim.y, .(.n, school), function(x) {
  m <- lm(y ~ standLRT, data = x)
  data.frame(x, fitted = fitted(m))
})

exam.true.df <- subset(Exam, select = c(school, normexam, standLRT))
colnames(exam.true.df)[2] <- "y"
exam.true.fitted <- ddply(exam.true.df, .(school), function(x) {
  m <- lm(y ~ standLRT, data = x)
  data.frame(x, fitted = fitted(m))	
})

# Creating the lineup: linear smoothers are plotted for each group,
# the different panes represent different data sets.

#      *********** Figure 2 ********************
qplot(x = standLRT, y = fitted, data = exam.true.fitted, group = school, 
      geom = "line", alpha = I(0.5)) %+% 
  lineup(true = exam.true.fitted, samples = exam.m1.fitted) + 
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 2 ********************


## Adding the random slope and refitting the LME model
exam.mod2 <- lmer(normexam ~ standLRT + (standLRT  | school), data = Exam)

# Simulating data for the null plots from model M1
exam.mod2.sims  <- simulate(exam.mod2, nsim = 19, seed = 1234)
exam.mod2.refit <- lapply(exam.mod2.sims, refit, object = exam.mod2)
exam.mod2.simy <- lapply(exam.mod2.refit, function(x) getME(exam.mod2, "y"))

exam.sim2.y <- do.call("cbind", exam.mod2.simy)
exam.sim2.y <- melt(exam.sim2.y)[,-1]
names(exam.sim2.y) <- c(".n", "y")
exam.sim2.y$.n <- as.numeric(str_extract(exam.sim2.y$.n, "\\d+"))
exam.sim2.y$standLRT <- rep(Exam$standLRT, rep = 19)
exam.sim2.y$school <- rep(Exam$school, rep = 19)

# Extracting information for the linear smoothers than will be plotted for each group
exam.mod2.fitted <- ddply(exam.sim2.y, .(.n, school), function(x) {
  m <- lm(y ~ standLRT, data = x)
  data.frame(x, fitted = fitted(m))
})

exam.true2.df <- subset(Exam, select = c(school, normexam, standLRT))
colnames(exam.true2.df)[2] <- "y"

exam.true2.fitted <- ddply(exam.true2.df, .(school), function(x) {
  m <- lm(y ~ standLRT, data = x)
  data.frame(x, fitted = fitted(m))	
})

#      *********** Figure 4 ********************
qplot(x = standLRT, y = fitted, data = exam.true2.fitted, group = school, 
      geom = "line", alpha = I(0.5)) %+% 
  lineup(true = exam.true2.fitted, samples = exam.mod2.fitted) + 
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 4 ********************


# Fitting a model with independent random effects
exam.mod3 <- lmer(normexam ~ standLRT + (standLRT - 1 | school) + (1 | school), data = Exam)

# Simulating random effects from model exam.mod3 with indepenent random effects
set.seed(987654321)
exam.mod3.sims  <- simulate(exam.mod3, nsim = 19)
exam.mod3.refit <- lapply(exam.mod3.sims, refit, object = exam.mod3)
exam.mod3.sim.ranef <- lapply(exam.mod3.refit, function(x) ranef(x)[[1]])

# Formatting a data frame of the null plots for plotting
exam.mod3.sim.ranef <- do.call("rbind", exam.mod3.sim.ranef)
exam.mod3.sim.ranef$.n <- rownames(exam.mod3.sim.ranef)
exam.mod3.sim.ranef$.n <- as.numeric(str_extract(exam.mod3.sim.ranef$.n, "\\d+"))

# Extracting the random effects from the model where we allow random 
# effects to be correlated
true.exam.mod2.ranef <- ranef(exam.mod2)$school

# Adding the true data to the data frame for the null plots
true_pos <- sample(20, 1) # print to reveal position
fig5_df <- nullabor:::add_true(exam.mod3.sim.ranef, true.exam.mod2.ranef, pos = true_pos)


#      *********** Figure 5 ********************
ggplot(data = fig5_df, aes(x = `(Intercept)`, y = standLRT)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.4) +
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 5 ********************

#-------------------------------------------------------------------------------
# Lineups based on the autism data set
# Figures 3 and 13
#-------------------------------------------------------------------------------

# The null model
autism.mod0 <- lmer(vsae ~ poly(age2, 2) + (age2 -1 | childid) + (I(age2^2) - 1 | childid), 
           data = autism)

# Simulating the level-1 residuals for the null plots
set.seed(9221632)
autism.mod0.sim.sicdegp  <- data.lineup.explvar1(null.model = autism.mod0, variable = "sicdegp", 
                                        data = autism, std = FALSE)

# Extracting the level-1 residuals for the true plot
autism.mod0.true.sicdegp <- data.frame(residual = resid(autism.mod0), sicdegp = autism$sicdegp)


#      *********** Figure 3 ********************
qplot(x = sicdegp, y = residual, data = autism.mod0.true.sicdegp, geom = "boxplot", 
      fill = sicdegp, outlier.size = 2, alpha=I(0.6)) %+% 
  lineup(true = autism.mod0.true.sicdegp, samples = autism.mod0.sim.sicdegp) + 
  facet_wrap( ~ .sample, ncol=5) + 
  ylim(-10, 10) + 
  ylab(NULL) + 
  xlab(NULL) + 
  scale_fill_brewer("", palette="Set2", labels=c("low", "medium", "high")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 3 ********************


# Model with only a linear random slope
autism.mod1 <- lmer(vsae ~ poly(age2, 2) + (age2 - 1 | childid), data = autism)

# Model with a quadratic random slope
autism.mod2 <- lmer(vsae ~ poly(age2, 2) + (age2 -1 | childid) + (I(age2^2) - 1 | childid), 
           data = autism)


# Simulating data for the null plots from autism.mod1
autism.mod1.sims  <- simulate(autism.mod1, nsim = 19, seed = 12345)
autism.mod1.refit <- lapply(autism.mod1.sims, refit, object = autism.mod1)
autism.mod1.sim.y <- lapply(autism.mod1.refit, function(x) getME(x, "y"))

# Formatting the simulated data for use in a lineup
autism.mod1.sim.y <- do.call("cbind", autism.mod1.sim.y)
autism.mod1.sim.y <- melt(autism.mod1.sim.y)[,-1]
names(autism.mod1.sim.y) <- c(".n", "y")
autism.mod1.sim.y$y[autism.mod1.sim.y$y < 0] <- 0
autism.mod1.sim.y$childid <- rep(autism$childid, rep = 19)
autism.mod1.sim.y$age2 <- rep(autism$age2, rep = 19)
autism.mod1.sim.y$.n <- as.numeric(str_extract(autism.mod1.sim.y$.n, "\\d+"))

# Formatting the true data
autism.true.y <- data.frame(y = autism$vsae, age2 = autism$age2, childid = autism$childid)

# Adding the true data to the data frame for the null plots
fig13_true_pos <- sample(20, 1) # print to reveal position
fig13_df <- nullabor:::add_true(autism.mod1.sim.y, autism.true.y, pos = fig13_true_pos)

#      *********** Figure 13 ********************
ggplot(data = fig13_df, aes(x = age2, y = y, group = childid)) +
  geom_line(alpha = 0.4) +
  facet_wrap( ~ .sample, ncol=5) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 13 ********************


#-------------------------------------------------------------------------------
# Lineups based on the Dialyzer data set
# Figures 6, 9, 12, 14 and 15
#-------------------------------------------------------------------------------

data("Dialyzer", package="MEMSS")

# Fitting the quartic LME model
dialyzer.mod2 <- lmer(rate ~ (pressure + I(pressure^2) + I(pressure^3) + I(pressure^4))*QB + 
             (pressure + I(pressure^2) | Subject), data = Dialyzer)

# Extracting the residuals
dialyzer.mod2.resid.df <- data.frame(dialyzer.mod2@frame, resid = resid(dialyzer.mod2))

# Simulating residuals for the null plots from model dialyzer.mod2
dialyzer.mod2.sims <- simulate(dialyzer.mod2, nsim = 19)
dialyzer.mod2.refit <- lapply(dialyzer.mod2.sims, refit, object = dialyzer.mod2)
dialyzer.mod2.sim.resids <- ldply(dialyzer.mod2.refit, function(x) data.frame(x@frame, resid = resid(x)))
dialyzer.mod2.sim.resids$.n <- as.numeric(str_extract(dialyzer.mod2.sim.resids$.id, "\\d+"))


#      *********** Figure 6 ********************
qplot(pressure, resid, data = dialyzer.mod2.resid.df,
      geom = "point", alpha = I(0.5)) %+%
  lineup(true = dialyzer.mod2.resid.df, samples = dialyzer.mod2.sim.resids) +
  facet_wrap(~ .sample, ncol = 5) +
  xlab(NULL) + 
  ylab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 6 ********************



# Fitting the quadratic model
dialyzer.mod1 <- lmer(rate ~ (pressure + I(pressure^2))*QB + (pressure + I(pressure^2) | Subject), 
           data = Dialyzer)

# Extracting the residuals
dialyzer.mod1.resid.df <- data.frame(dialyzer.mod1@frame, resid = resid(dialyzer.mod1))

# Simulating residuals for the null plots from model dialyzer.mod1
set.seed("20040110")
dialyzer.mod1.sims <- simulate(dialyzer.mod1, nsim = 19)
dialyzer.mod1.refit <- lapply(dialyzer.mod1.sims, refit, object = dialyzer.mod1)
dialyzer.mod1.sim.resids <- ldply(dialyzer.mod1.refit, function(x) data.frame(x@frame, resid = resid(x)))
dialyzer.mod1.sim.resids$.n <- as.numeric(str_extract(dialyzer.mod1.sim.resids$.id, "\\d+"))

dialyzer.mod1.resids <- rbind.fill(dialyzer.mod1.resid.df, dialyzer.mod1.sim.resids)

# Randomizing the panels
fig9.plot.order <- sample.int(20, 20)
fig9.location <- fig9.plot.order[1] # location of the true plot
dialyzer.mod1.resids$.id <- rep(fig9.plot.order, each = nrow(dialyzer.mod1@frame))


#      *********** Figure 9 ********************
ggplot(dialyzer.mod1.resids, aes(x = pressure, y = resid)) + 
  geom_point() + 
  geom_smooth(method = "loess") +
  facet_wrap(~ .id, ncol = 5) +
  xlab(NULL) + 
  ylab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 9 ********************


# Discretizing pressure for boxplots representation
fig12_df <- dialyzer.mod2.sim.resids
fig12_df$pressurecat <- cut(fig12_df$pressure, 
                            breaks = c(0, 0.45, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25))

fig12_true <- dialyzer.mod2.resid.df
fig12_true$pressurecat <- cut(fig12_true$pressure, breaks = c(0, 0.45, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25))

# Adding the true data to the data frame for the null plots
fig12_true_pos <- sample(20, 1) # print to reveal position
fig12_df <- nullabor:::add_true(fig12_df, fig12_true, pos = true_pos)

#      *********** Figure 12 ********************
qplot(x = pressurecat, y = resid, data = fig12_df, geom = "boxplot", fill = pressurecat, alpha=I(0.8)) +
  facet_wrap(~ .sample, ncol = 5) +
  xlab(NULL) + 
  ylab(NULL) + 
  scale_fill_brewer() +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.position="none") + 
  scale_x_discrete(expand = c(.1,.1))
#      *********** Figure 12 ********************



# Formatting data frame for plotting and sorting Subjects by residual variance
fig14_true_pos <- sample(20, 1) # print to reveal position
fig14_df <- nullabor:::add_true(dialyzer.mod1.sim.resids, dialyzer.mod1.resid.df, pos = fig14_true_pos)

fig14_df <- ddply(fig14_df, .(.sample, Subject), transform, var = var(resid))
fig14_df <- ddply(fig14_df, .(.sample), transform, subid = rank(var, ties.method="min"))

#      *********** Figure 14 ********************
qplot(subid, resid, data = fig14_df, geom = c("point"), fill=I("grey70"), alpha=I(0.5)) +
  facet_wrap(~ .sample, ncol = 5) +
  geom_hline(yintercept=0, colour="red") +
  xlab(NULL) + 
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.position="none") +
  scale_x_continuous(expand = c(.075,.075))
#      *********** Figure 14 ********************


#      *********** Figure 15 ********************
qplot(Subject, resid, data = dialyzer.mod1.resid.df, geom = "boxplot", fill=I("grey70"), alpha=I(0.8)) %+%
  lineup(true = dialyzer.mod1.resid.df, samples = dialyzer.mod1.sim.resids) +
  facet_wrap(~ .sample, ncol = 5) +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position="none") +
  scale_x_discrete(expand = c(.075,.075))
#      *********** Figure 15 ********************


#-------------------------------------------------------------------------------
# Lineup based on the ahd data
# Figure 7
#-------------------------------------------------------------------------------

# Fitting the proposed model
ahd.model <- lmer(sbvalue ~ treatment + week*baseline + I(week^2)*baseline + 
                (0 + week + I(week^2) | subject), data = ahd)

# Checking between group homogeneity
## Residuals from true model
set.seed(123456789)

## Simulating null data and refitting the model
ahd.nsim <- 19
ahd_refit_fm <- vector(mode = "list", length = ahd.nsim)
for(i in 1:ahd.nsim) {
  ahd_refit_fm[[i]] <- refit(ahd.model, simulate(ahd.model))
}

## Extract level-1 residuals both from the sims and the true model
ahd_sim_resids <- ldply(ahd_refit_fm, function(x){
  e <- resid(x)
  subj <- x@frame$subject
  return(data.frame(EB.resid = e, subject = subj))
})
ahd_sim_resids <- data.frame(.id = rep(1:19, each = nrow(ahd.model@frame)), ahd_sim_resids)

ahd_true_resids <- data.frame(.id = 20, EB.resid = resid(ahd.model), subject = ahd.model@frame$subject)

# combining the results for plotting
ahd_resids <- rbind(ahd_sim_resids, ahd_true_resids)

# randomizing the order of the plots within the lineup
ahd_resids$sample <- sample(20, 20, replace = FALSE)[ahd_resids$.id]

# Ordering the subjects by IQR
ahd_resids <- ddply(ahd_resids, .(.id, subject), transform, iqr = IQR(EB.resid))
ahd_resids <- ddply(ahd_resids, .(.id), transform, rank = order(order(iqr, subject)))
ahd_resids <- ddply(ahd_resids, .(.id, subject), transform, rank = min(rank))
ahd_resids <- ddply(ahd_resids, .(.id), transform, rank = rank(rank))

#      *********** Figure 7 ********************
qplot(x = factor(rank), y = EB.resid, data = ahd_resids, 
      geom = "boxplot", xlab = "", ylab = "", outlier.size = 1.5) + 
  coord_flip() + 
  ylim(-125, 125) + 
  theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank(), panel.grid.major.y = element_blank()) +
  facet_wrap(~sample)
#      *********** Figure 7 ********************



#-------------------------------------------------------------------------------
# Lineups for the radon data
# Figures 8, 10
#-------------------------------------------------------------------------------

# Calculating the no. of observations for each county
radon <- ddply(radon, .(county), transform, n = length(county))

# The model using only counties with at least 4 obs.
radon_fm <- lmer(log.radon ~ basement + uranium + (basement | county), 
           data = subset(radon, n > 4))


# The *observed residuals* and other variables needed for lineup
observed.cyclone.df <- data.frame(radon_fm@frame, resid = resid(radon_fm))

# Simulating a set of 19 null plots
radon.sim.y   <- simulate(radon_fm, nsim = 19, seed = 987654321) ## Change seed for diff. nulls                       
radon.sim.mod <- apply(radon.sim.y, 2, refit, object = radon_fm)

# Extract level-1 residuals both from the simulated models
null.cyclone.df <- ldply(radon.sim.mod, function(x) data.frame(x@frame, resid = resid(x)))
null.cyclone.df$.id <- as.numeric( str_extract(null.cyclone.df$.id, "\\d+") )

# Combine for plotting
cyclone.df <- rbind.fill(observed.cyclone.df, null.cyclone.df)

# Randomizing the panels
fig8.plot.order <- sample.int(20, 20)
fig8.location <- fig8.plot.order[1] # location of the true plot
cyclone.df$.id <- rep(fig8.plot.order, each = nrow(radon_fm@frame))

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


#      *********** Figure 8 ********************
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
#      *********** Figure 8 ********************


## Q-Q plots

# Fit the proposed model
radon_fm2 <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

# Extract the random slope
radon_fm2_b <- ranef(radon_fm2)[[1]] # notice that this is actually a matrix

# Simulating data from the model
radon2.sim.y   <- simulate(radon_fm2, nsim = 19, seed = 987654321)                        
radon2.sim.mod <- apply(radon2.sim.y, 2, refit, object = radon_fm2)            ## a list of models

# Extracting the random slopes from the simulations
# and formating a data set of the 20 sets of Q-Q plots
radon2.sim.b1 <- llply(radon2.sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
radon2.sim.b1 <- melt( do.call("rbind", radon2.sim.b1) )[,-2]           ## changing to a data frame
names(radon2.sim.b1) <- c("sample", "basement")                         ## setting colnames for faceting
radon2.sim.b1        <- arrange(radon2.sim.b1, sample)                  ## ordering by simulation

radon2.sim.b1$.n <- as.numeric( str_extract(radon2.sim.b1$sample, "\\d+") )
radon2.sim.b1 <- ddply(radon2.sim.b1[complete.cases(radon2.sim.b1),], .(.n), transform, 
                       band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

radon2.b1 <- transform(radon_fm2_b, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

#      *********** Figure 10 ********************
qplot(sample = basement, data = radon2.b1, stat = "qq") %+%
  lineup(true = radon2.b1, sample = radon2.sim.b1) + 
  facet_wrap(~ .sample, ncol = 5) + 
  geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
  ylab(NULL) + 
  xlab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 10 ********************

#-------------------------------------------------------------------------------
# Figure 11: Lineup of Q-Q plots for the random slope simulated from a t dsn
#-------------------------------------------------------------------------------

# Fitting the proposed model
fig11_fm <- lmer(log.radon ~ basement + uranium + (1 | county) + 
             (basement - 1 | county), data = radon)

# Simulating data with random effects from a t distribution
set.seed(-2029298609)
fig11.y.b.t <- sim_t_hlm(fig11_fm)
fig11.refit.b.t <-  refit(fig11_fm, fig11.y.b.t)
fig11.b.t <- ranef(fig11.refit.b.t)[[1]]

# Simulating a normal LME model
fig11.sim.y   <- simulate(fig11_fm, nsim = 19)                        
fig11.sim.mod <- apply(fig11.sim.y, 2, refit, object = fig11_fm)  ## a list of models

# Extracting the random slopes from the simulations
# and formating a data set of the 20 sets of Q-Q plots
fig11.sim.b1 <- llply(fig11.sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
fig11.sim.b1 <- melt( do.call("rbind", fig11.sim.b1) )[,-2]           ## changing to a data frame
names(fig11.sim.b1) <- c("sample", "basement")                        ## setting colnames for faceting
fig11.sim.b1        <- arrange(fig11.sim.b1, sample)                  ## ordering by simulation

fig11.sim.b1$.n <- as.numeric( str_extract(fig11.sim.b1$sample, "\\d+") )
fig11.sim.b1 <- ddply(fig11.sim.b1, .(.n), transform, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

fig11.b1 <- transform(fig11.b.t, band = sim_env(basement), 
                x = sort(qqnorm(basement, plot.it=FALSE)$x))

#      *********** Figure 11 ********************
fig11.location <- 6
qplot(sample = basement, data = fig11.b1, stat = "qq") %+%
  lineup(true = fig11.b1, sample = fig11.sim.b1, pos=fig11.location) + 
  facet_wrap(~ .sample, ncol = 5) + 
  geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
  #	xlab("Normal Quantiles") + ylab("Sample Quantiles") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
#      *********** Figure 11 ********************


#-------------------------------------------------------------------------------
# Table 1
#-------------------------------------------------------------------------------

# Exploring the conventional test for homoscedastic error 
# terms across groups for different minimum group sizes.

# Calculating no. of obs. per county
radon <- ddply(radon, .(county), transform, n = length(county))

# Using different cutoffs for what group size is considered "small"
RES <- NULL
for(i in 3:15) {
  RES <- rbind(RES, c(n = i, het.chisq.test(cutoff = i)))
}
RES <- as.data.frame(RES)


# A simulation-based version of the conventional test for homoscedastic error 
# terms across groups using the parametric bootstrap. Different minimum group 
# sizes are considered.

RES$sim.p.val <- rep(NA, 13)
for(i in 1:13) {
  set.seed(987654321)
  fm <- lmer(log.radon ~ basement + uranium + (basement | county), 
             data = subset(radon, n >= RES$n[i]))
  sim.ys <- simulate(fm, nsim = 1e4)
  sim.df <- lapply(sim.ys, function(y) {
    df <- fm@frame
    df$log.radon <- y
    return(df)
  })
  
  sim.h <- mclapply(sim.df, FUN = function(x){
    df <- ddply(x,  .(county), calc_s_df, formula = log.radon ~ basement)
    df <- transform(df, d = ( log(s^2) - ( sum(df * log(s^2)) / sum(df) ) ) / sqrt( 2 / df ) )
    return(sum(df$d^2))
  })
  
  RES$sim.p.val[i] <- mean(sim.h >= RES$H[i])
  cat("Iteration", i, "complete \n", sep = " ")
}

#      *********** Table 1 ********************
print(xtable(RES), sanitize.text.function=function(x)x)
#      *********** Table 1 ********************



#-------------------------------------------------------------------------------
# Creating the tables from the Supplement
#-------------------------------------------------------------------------------

turk <- read.csv("../data/study.csv")
turk <- turk %>% dplyr::mutate(
  correct = response==data_location
)


lps <- ddply(turk, .(lineup), summarise,
             correct=sum(correct*weight),
             num=sum(weight)
)

lps$pvals <- unlist(llply(1:nrow(lps), function(i) {
  correct <- floor(lps$correct[i])
  n <- lps$num[i]
  vinference:::pV(correct, n, m=20, scenario=3)
}))
lps$signif <- lps$pvals < 0.05

lps$rep <- as.numeric(gsub(".*-([0-9])","\\1",as.character(lps$lineup)))
lps$exp <- gsub("(.*)-[0-9]","\\1",as.character(lps$lineup))

pval.overall <- ddply(lps, .(exp), function(x) {
  res <- vinference:::scenario3(N=10000, sum(x$num))
  dres <- as.data.frame(res)
  with(dres, sum(Freq[as.numeric(names(res)) >= sum(x$correct)]))
})$V1



lps$stars <- cut(lps$pvals, breaks=c(0,0.001, 0.01, 0.05, 0.1, 1))
levels(lps$stars) <- c("***", "**", "*", ".", " ")

lps$str <- with(lps, sprintf("%d/%d & \\hspace{-0.1in}%s", floor(correct), round(num), stars))

dt <- dcast(lps, exp~rep, value.var="str")

dt$overall <- ifelse(pval.overall < 10^-4, "$< 10^{-4}$", sprintf("%.4f",pval.overall))


#      *********** Table 1 Supplement ********************
print(xtable(dt), sanitize.text.function=function(x)x)
#      *********** Table 1 Supplement ********************


turk$choice <- gsub("_.*", "", as.character(turk$reason))
turk$choiceWT <- nchar(turk$choice)

reasons <- dlply(turk, .(lineup, correct), function(x) {
  choices <- unlist(strsplit(x$choice, split=""))
  weights <- rep(x$weight*1/x$choiceWT, x$choiceWT)
  dt <- xtabs(weights~choices)
  as.data.frame(dt)
})
dreasons <- ldply(reasons, function(x) x)
dreasons$exp <- gsub("(.*)-[0-9]", "\\1", as.character(dreasons$lineup))
dreasons$pick <- c("null", "data")[dreasons$correct+1]

res <- ddply(dreasons, .(exp, pick, choices), summarise, Freq=sum(Freq))


# probabilities the other way round
qt <- ddply(res, .(exp, choices), transform, perc=Freq/sum(Freq)*100)
qt2 <- dcast(qt, exp+pick~choices, value.var="perc")
names(qt2)[3:7] <- c("Outlier", "Spread", "Trend", "Asymmetry", "Other")

#      *********** Table 2 Supplement ********************
print(xtable(subset(qt2, pick=="data")[,-2], digits=c(1,1,1,1,1,1,1)), include.rownames=FALSE, NA.string="0.0")
#      *********** Table 2 Supplement ********************