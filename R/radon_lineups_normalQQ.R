#-------------------------------------------------------------------------------
# Script fitting the model(s) to the radon data; and obtaining residuals
# for lineups. 
#
# Adam Loy
# August 2013
#-------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
### Preliminaries
# -----------------------------------------------------------------------------------------
library(lme4)
library(plyr)
library(reshape2)
library(stringr)

load(file.choose()) ## find radon_for_sims.csv

### Function to reorganize Z like it should be in 511
BlockZ <- function(object) {
  Z <- getME(object, "Z")
  
  grp.size <- table(object@flist)
  ngrps <- length(grp.size)
  nranef <- dim(ranef(object)[[1]])[2]
  
  base.ord <- seq(from = 1, by = ngrps, length.out = nranef)
  ord <- base.ord + rep(0:(ngrps - 1), each = nranef)
  
  perm.mat <- t(as(ord, "pMatrix"))
  
  return(Z %*% perm.mat)
}

### Function to calculate the marginal variance of random effects
lev2.marginal.var <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  ngrps <- unname(sapply(.model@flist, function(x) length(levels(x))))
  
  # Constructing V = Cov(Y)
  sig0 <- attr(VarCorr(.model), "sc") # sigma(.model)
  
  ZDZt <- sig0^2 * crossprod( .model@A )
  R    <- Diagonal( n = n, x = sig0^2 )
  D    <- kronecker( Diagonal(ngrps), bdiag(VarCorr(.model)) )
  V    <- Diagonal(n) + ZDZt
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )

  bse <- crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
  bse.diag <- diag(bse)

  semat <- matrix(sqrt(bse.diag), ncol = 2, byrow = TRUE)

  return(semat)
}

### Function to extract standardized random effects
std_ranef <- function(.model) {
	res <- ranef(.model)[[1]]
	semat <- lev2.marginal.var(.model)
	
	RVAL <- res / semat
	return(RVAL)
}

### Function for weighted empirical CDF
wecdf <- function(x, weights) {
    stopifnot(length(x) == length(weights))
    sw <- sum(weights)
    if (length(x) < 1) 
        stop("'x' must have 1 or more non-missing values")
    stopifnot(all(weights >= 0))
    ox <- order(x)
    x  <- x[ox]
    w  <- weights[ox]
    vals <- sort(unique(x))
    xmatch <- factor(match(x, vals), levels = seq_along(vals))
    wmatch <- tapply(w, xmatch, sum)
    wmatch[is.na(wmatch)] <- 0
    rval <- approxfun(vals, cumsum(wmatch) / sw, method = "constant", 
        yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    return(rval)
}   



# -----------------------------------------------------------------------------------------
### Fitting the models
# -----------------------------------------------------------------------------------------

# Random intercept, random slope model
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = radon)

# -----------------------------------------------------------------------------------------
### Calculating the original residuals
# -----------------------------------------------------------------------------------------

# Level-1 residuals (i.e., the error terms)
e <- resid(fm)

# Random effects (i.e., the level-2 residuals)
b <- ranef(fm)[[1]] # notice that this is actually a matrix

# -----------------------------------------------------------------------------------------
### Simulating null residuals via the parametric bootstrap
# -----------------------------------------------------------------------------------------
set.seed(9221632)

for(i in 1:5) {

sim.y   <- simulate(fm, nsim = 19)                        ## A 919 x 19 matrix of responses
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models

# Simulated random slopes
sim.b1 <- ldply(sim.mod, function(x) {                    ## a df of random slopes
	data.frame(slope = ranef(x)[[1]][,2], 
	           std.slope = std_ranef(x)[,2])
}) 
names(sim.b1)[1] <- "sample"                              ## setting colnames for faceting
# sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation
sim.b1.df <- rbind(sim.b1,
			 data.frame(sample = "true", slope = b[,2],
			       std.slope = std_ranef(fm)[,2]))        ## adding the true r. slopes

# saving the 5 replications to different objects
nam <- paste("sim_slope_qq_df", i, sep = "_")
assign(nam, sim.b1.df)
}

# I named the file: 
# save(sim_slope_qq_df_1, sim_slope_qq_df_2, sim_slope_qq_df_3, 
#      sim_slope_qq_df_4, sim_slope_qq_df_5, file = file.choose(new = T))