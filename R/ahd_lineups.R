#-------------------------------------------------------------------------------
# Script constructing lineups of the cyclone plots.
#
# Adam Loy
# May 2013
#
# METHYLPREDISONE DATA
# Description: 
# 
#
# References: - Carithers et al. (1989). Methylprednisolone therapy in patients
#               with severe alcoholic hepatitis. Annals of Internal Medicine, 
#               110(9):685–690.
#             - Vonesh, E. F and Chinchilli, V. M. (1997). Linear and Nonlinear
#               Models for the Analysis of Repeated Measurements. Marcel 
#               Dekker, New York.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------
setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")
load('data/ahd.RData') # data
library(lme4)
library(HLMdiag)
library(ggplot2)
library(grid)
library(plyr)
library(MASS)

model <- lmer(sbvalue ~ treatment + week*baseline + I(week^2)*baseline + 
  (0 + week + I(week^2) | subject), data = ahd)


#-------------------------------------------------------------------------------
# Cyclone plots - Assumption VIOLATED
#-------------------------------------------------------------------------------

# Checking between group homogeneity
## Residuals from true model
set.seed(123456789)

## Simulating null data
sims <- simulate(model, nsim = 19)

## Refit
refit_fm <- apply(sims, 2, refit, object = model)

## Extract level-1 residuals
sim_resids <- llply(refit_fm, function(x){
	e <- resid(x)
	subj <- x@frame$subject
	return(data.frame(EB.resid = e, subject = subj))
})

true_resids <- data.frame(EB.resid = resid(model), subject = model@frame$subject)

pdf.options(reset = FALSE)
pdf("figures/ahd_badcyclone5.pdf", width = 8.5, height = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5,4)))
vplayout <- function(x, y){viewport(layout.pos.row = x, layout.pos.col = y)}
realp <- qplot(x = reorder(subject, EB.resid, IQR), y = EB.resid, data = true_resids, 
               geom = "boxplot", xlab = "subject", ylab = "residuals") + coord_flip() + 
                 ylim(-150, 150) + 
                 theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank())
real.i <- 5 # sample(1:20, 1)
j <- 0
pos <- matrix(1:20, ncol = 4, byrow=T)
for(i in 1:20){
  if(i==real.i) { 
    if(!i %in% c(1,5,9,13,17)) realp <- realp + xlab(NULL)
    if(!i %in% 17:20) realp <- realp + ylab(NULL)
    
    print(realp, vp = vplayout(which(pos == i, arr.ind = TRUE)[1],
                               which(pos == i, arr.ind = TRUE)[2]))
  }
  else{
    j <- j + 1
    p <- qplot(x = reorder(subject, EB.resid, IQR), y = EB.resid, 
               data = sim_resids[[j]], geom = "boxplot", xlab = "subject", ylab = "residuals") + 
                 coord_flip() + ylim(-150, 150) + 
                 theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank())
    
    if(!i %in% c(1,5,9,13,17)) p <- p + xlab(NULL)
    if(!i %in% 17:20) p <- p + ylab(NULL)
    
    print(p, vp = vplayout(which(pos == i, arr.ind = TRUE)[1], which(pos == i, arr.ind = TRUE)[2]))
  }
}
dev.off()

#-------------------------------------------------------------------------------
# Cyclone plots - Assumption OK
#-------------------------------------------------------------------------------

model2 <- lmer(log(sbvalue) ~ treatment + week*baseline + I(week^2)*baseline + 
  (0 + week + I(week^2) | subject), data = ahd)


## Simulating null data
sims2 <- simulate(model2, nsim = 19)

## Refit
refit_fm2 <- apply(sims2, 2, refit, object = model2)

## Extract level-1 residuals
sim_resids2 <- llply(refit_fm2, function(x){
	e <- resid(x)
	subj <- x@frame$subject
	return(data.frame(EB.resid = e, subject = subj))
})

true_resids2 <- data.frame(EB.resid = resid(model2), subject = model2@frame$subject)

pdf.options(reset = FALSE)
pdf("figures/ahd_goodcyclone13.pdf", width = 8.5, height = 11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5,4)))
vplayout <- function(x, y){viewport(layout.pos.row = x, layout.pos.col = y)}
realp <- qplot(x = reorder(subject, EB.resid, IQR), y = EB.resid, data = true_resids2, 
               geom = "boxplot", xlab = "subject", ylab = "residuals") + coord_flip() + 
                 ylim(-.5, .5) + 
                 theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank())
real.i <- 5 # sample(1:20, 1)
j <- 0
pos <- matrix(1:20, ncol = 4, byrow=T)
for(i in 1:20){
  if(i==real.i) { 
    if(!i %in% c(1,5,9,13,17)) realp <- realp + xlab(NULL)
    if(!i %in% 17:20) realp <- realp + ylab(NULL)
    
    print(realp, vp = vplayout(which(pos == i, arr.ind = TRUE)[1],
                               which(pos == i, arr.ind = TRUE)[2]))
  }
  else{
    j <- j + 1
    p <- qplot(x = reorder(subject, EB.resid, IQR), y = EB.resid, 
               data = sim_resids2[[j]], geom = "boxplot", xlab = "subject", ylab = "residuals") + 
                 coord_flip() + ylim(-.5, .5) + 
                 theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank())
    
    if(!i %in% c(1,5,9,13,17)) p <- p + xlab(NULL)
    if(!i %in% 17:20) p <- p + ylab(NULL)
    
    print(p, vp = vplayout(which(pos == i, arr.ind = TRUE)[1], which(pos == i, arr.ind = TRUE)[2]))
  }
}
dev.off()