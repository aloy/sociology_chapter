### Performing the "classical" analogs to the lineups
### August 2014

library(lme4)
library(plyr)

### Selection of random effects
### Data set: Exam
data("Exam", package = "mlmRev")

# Likelihood ratio test for a random slope
M0 <- lmer(normexam ~ standLRT + (1  | school), data = Exam)
M1 <- lmer(normexam ~ standLRT + (standLRT  | school), data = Exam)
anova(M1, M0) # comparing to chisq w/ 2 df
mean(pchisq(40.372, df = 1:2, lower.tail = FALSE)) # using 50:50 mixture

# Likelihood ratio test for the correlation between random effects
M2 <- lmer(normexam ~ standLRT + (1  | school) + (0 + standLRT | school), data = Exam)
anova(M1, M2)



### Checking homogeneity between groups
### Data set: ahd
data("ahd", package = "HLMdiag")

# Fitted model
fm <- lmer(sbvalue ~ treatment + week*baseline + I(week^2)*baseline + 
             (0 + week + I(week^2) | subject), data = ahd)


# Test suggested by Raudenbush and Byrk
# To perform the test we need to eliminate subjects with less than 4 observations,
# though it is quite dubious to base anything off models fit to such little data!
nobs.subj <- ddply(ahd, .(subject), summarise, nobs = sum(!is.na(sbvalue)))
few.obs <- nobs.subj[nobs.subj$nobs < 4, 1]
reduced.ahd <- subset(ahd, !subject %in% few.obs)
  
# sep.lms <- lmList(sbvalue ~ week + I(week^2) | subject, data = reduced.ahd)
calc_s_df <- function(x, formula) {
  mod <- lm(formula, data = x)
  smry <- summary(mod)
  s  <- smry$sigma  # residual std. dev.
  df <- smry$df[2]  # n_i - r_i
  
  return(c(s = s, df = df))
  #   res <- ( log(s^2) - ( sum(df * log(s^2)) / sum(df) ) ) / sqrt( 2 / df )
  #   return(res)
}

test.df <- ddply(reduced.ahd, .(subject), calc_s_df, formula = sbvalue ~ week + I(week^2))
test.df <- transform(test.df, d = ( log(s^2) - ( sum(df * log(s^2)) / sum(df) ) ) / sqrt( 2 / df ) )

H <- sum(test.df$d^2)
pchisq(H, df = nrow(test.df) - 1, lower.tail = FALSE)

# Making sure I coded the above correctly
# ls.pooled <- (1/sum(test.df$df)) * sum(test.df$df * log(test.df$s^2))
# H <- sum( (test.df$df/2) * (log(test.df$s^2) - ls.pooled)^2 )



### Checking homogeneity between groups
### Data set: radon
data("radon", package = "HLMdiag")

radon <- ddply(radon, .(county), transform, n = length(county))

# Model from which cyclone lineup was created
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = subset(radon, n > 9))

test.df2 <- ddply(fm@frame, .(county), calc_s_df, formula = log.radon ~ basement)
test.df2 <- transform(test.df, d = ( log(s^2) - ( sum(df * log(s^2)) / sum(df) ) ) / sqrt( 2 / df ) )

H2 <- sum(test.df2$d^2)
pchisq(H2, df = nrow(test.df) - 1, lower.tail = FALSE)


# Checking the p-value via simulation

set.seed(987654321)
sim.ys <- simulate(fm, nsim = 1000)
sim.df <- lapply(sim.ys, function(y) {
  df <- fm@frame
  df$log.radon <- y
  return(df)
  })

sim.h <- sapply(sim.df, FUN = function(x){
  df <- ddply(x,  .(county), calc_s_df, formula = log.radon ~ basement)
  df <- transform(df, d = ( log(s^2) - ( sum(df * log(s^2)) / sum(df) ) ) / sqrt( 2 / df ) )
  return(sum(df$d^2))
})