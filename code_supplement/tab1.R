#-------------------------------------------------------------------------------
# Script generating table 1
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

library(checkpoint)
checkpoint("2016-06-10")

library(lme4)
library(plyr)
library(parallel)


# Loading the radon data set
data("radon", package = "HLMdiag")

# Calculating no. of obs. per county
radon <- ddply(radon, .(county), transform, n = length(county))

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
# Exploring the conventional test for homoscedastic error 
# terms across groups for different minimum group sizes.
#-------------------------------------------------------------------------------

# Using different cutoffs for what group size is considered "small"
RES <- NULL
for(i in 3:15) {
  RES <- rbind(RES, c(n = i, het.chisq.test(cutoff = i)))
}
RES <- as.data.frame(RES)

#-------------------------------------------------------------------------------
# A simulation-based version of the conventional test for homoscedastic error 
# terms across groups using the parametric bootstrap. Different minimum group 
# sizes are considered.
#-------------------------------------------------------------------------------

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
  
  RES$sim.p.val[i] <- mean(sim.h >= RES$H)
  cat("Iteration", i, "complete \n", sep = " ")
}

# Print to see the results
RES