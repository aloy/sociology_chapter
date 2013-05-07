#-------------------------------------------------------------------------------
# Script constructing lineups for use in a paper demonstrating their utility
# in checking the assumptions for hierarchical linear models.
#
# Adam Loy
# April 2013
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")

library(lme4)
library(HLMdiag)
library(foreign)
library(stringr)
library(HLMdiag)
library(nullabor)
library(plyr)
library(reshape2)

#-------------------------------------------------------------------------------
# An example of the usefulness of lineups in model exploration/selection we
# can use PISA 2009 data from the United States, or even subset it further to the 
# Midwest. One such subset is made available by Snijders and Bosker (2012)
# http://www.stats.ox.ac.uk/~snijders/mlbook.htm
#
# Data description: 
# PISA 2009 - Contains data relating to students' reading proficiency. 
#
# Variables:
# *STRATUM - US region and Private vs. Public
# *ST01Q01 - Grade
# *ST04Q01 - Sex
# *AGE - Age
# *IMMIG - Immigration status (1 = native; 2 = second generation; 3 = first generation)
# *ESCS - Index of economic, social, and cultural status
# *STUDREL - Teacher student relations
# *JOYREAD - Joy/like reading
# *METASUM - Meta-cognition: summarizing
# *UNDREM - Meta-cognition: understanding and remembering
# *SC02Q01 - Pulic or Private school
# *SC04Q01 - School Community
# *SCHSIZE - Total school enrollment
# *SELSCH - Index of academic school selectivity
# *STRATIO - Student-teacher ratio
#-------------------------------------------------------------------------------

# function to sub NA for -999
makeNA <- function(x) ifelse(x == -999, NA, x)

# Reading in data
pisa <- read.spss("data/combineusa-999_c.sav", to.data.frame=TRUE)
head(pisa)

# insert NAs
pisa$IMMIG <- makeNA(pisa$IMMIG)
pisa$ESCS <- makeNA(pisa$ESCS)
pisa$STUDREL <- makeNA(pisa$STUDREL)
pisa$JOYREAD <- makeNA(pisa$JOYREAD)
pisa$METASUM <- makeNA(pisa$METASUM)
pisa$UNDREM <- makeNA(pisa$UNDREM)
pisa$SCHSIZE <- makeNA(pisa$SCHSIZE)

# Formatting some variable names
pisa$sex <- pisa$ST04Q01
pisa$grade <- pisa$ST01Q01

# Making a region variable
pisa$region <- sapply( str_split(pisa$STRATUM, pattern = " "), function(x) x[3])

# making subset for the midwest
midwest <- subset(pisa, subset = region == "Midwest")
midwest <- na.omit(midwest)

### Some initial exploration of useful plots

null <- lmer(METASUM ~ 1 + (1 | SCHOOLID), data = midwest)

( mod1 <- lmer(METASUM ~ sex + AGE + ESCS + IMMIG + (1 | SCHOOLID), data = midwest) )

qplot(x = midwest$ESCS, y = resid(mod1))
qplot(x = reorder(midwest$SCHOOLID, X=resid(mod1), FUN=IQR), y = resid(mod1), geom = "boxplot") + coord_flip()

#-------------------------------------------------------------------------------
# An example of the usefulness of lineups in model exploration/selection we
# can use PISA 2009 data from the United States, or even subset it further to the 
# Midwest. One such subset is made available by Snijders and Bosker (2012)
# http://www.stats.ox.ac.uk/~snijders/mlbook.htm
#
# This example looks at 8th grade students (~ 11 years old) in schools in the
# Netherlands. The variables include:
# *schoolnr: school id
# *pupilNR_new: pupil id
# *langPOST: score on language test
# *ses: socio-economic status
# *IQ_verb: verbal IQ
# *sex: sex of the student
# *Minority: 
# *denomina:
# *sch_ses:
# *sch_iqv:
# *sch_min:
#-------------------------------------------------------------------------------

mlbook_red <- read.table("data/mlbook2_r.dat", header=TRUE)
names(mlbook_red)

(M1 <- lmer(langPOST ~ IQ_verb + (1 | schoolnr), data = mlbook_red))

qplot(x = fitted(M1), y = resid(M1))
qplot(x = mlbook_red$IQ_verb, y = resid(M1))
qplot(x = mlbook_red$sch_iqv, y = resid(M1)) + geom_smooth(method="lm")

# Is ses useful?
(M2 <- lmer(langPOST ~ IQ_verb + sch_iqv + (1 | schoolnr), data = mlbook_red))
qplot(x = mlbook_red$ses, y = resid(M2)) + geom_smooth(method="lm")

# Do we need a random slope?
(M3 <- lmer(langPOST ~ IQ_verb + sch_iqv + ses + (IQ_verb | schoolnr), data = mlbook_red))

# Do we need to include gender?
(M4 <- lmer(langPOST ~ IQ_verb * ses + sch_iqv + (IQ_verb | schoolnr), data = mlbook_red))
qplot(x = IQ_verb * ses, y = resid(M3), data = mlbook_red) + geom_smooth(method="lm")

# Are the level-1 residual heteroscedastic?
qplot(x = reorder(mlbook_red$schoolnr, X=resid(M4), FUN=IQR), y = resid(M4), geom = "boxplot") + coord_flip()

# Are the level-1 residual heteroscedastic wrt IQ?
qplot(x = mlbook_red$IQ_verb, y = resid(M4))

#-------------------------------------------------------------------------------
# An example of the usefulness of lineups in model exploration/selection 
#
# We could use the exam data again that we explored with the JSS paper.
#
# Data description:
# 4059 students nested within 65 schools
# Variables - 
# *school - school id
# *student - student id
# *normexam - student's standardized exam score at 16
# *schgend - school's gender
# *schavg - 
# *vr -
# *intake -
# *standLRT - student's standardized score on the London reading test (age 11)
# *sex - student's gender
#-------------------------------------------------------------------------------
library(mlmRev)

# Initial model
M1 <- lmer(normexam ~ standLRT + (1 | school), data = Exam)

# Do we need a random slope?
qplot(x = standLRT, y = normexam, data = Exam, geom = "smooth", group = school, se = F, method = "lm")

m1.sims  <- simulate(M1, nsim = 19, seed = 1234)
m1.refit <- lapply(m1.sims, refit, object = M1)
m1.simy <- lapply(m1.refit, function(x) x@y)

sim.y <- do.call("cbind", m1.simy)
sim.y <- melt(sim.y)[,-1]
names(sim.y) <- c(".n", "y")
sim.y$.n <- as.numeric(str_extract(sim.y$.n, "\\d+"))
sim.y$standLRT <- rep(Exam$standLRT, rep = 19)
sim.y$school <- rep(Exam$school, rep = 19)

true.y <- data.frame(y = Exam$normexam, standLRT = Exam$standLRT, school = Exam$school)

qplot(x = standLRT, y = y, data = true.y, group = school, geom = "smooth", method = "lm", se=F) %+% lineup(true = true.y, samples = sim.y) + facet_wrap( ~ .sample, ncol=5)


# Including a random slope
M2 <- lmer(normexam ~ standLRT + (standLRT | school), data = Exam)

# Do other variables help?
qplot(x = standLRT, y = resid(M2), data = Exam, geom = c("point", "smooth"))

m2.sims  <- simulate(M2, nsim = 19, seed = 1234)
m2.refit <- lapply(m2.sims, refit, object = M2)
m2.resid <- lapply(m2.refit, resid)

sim.resid2 <- do.call("cbind", m2.resid)
sim.resid2 <- melt(sim.resid2)[,-1]
names(sim.resid2) <- c(".n", "resid")
sim.resid2$.n <- as.numeric(str_extract(sim.resid2$.n, "\\d+"))
sim.resid2$standLRT <- rep(Exam$standLRT, times=19)

true.resid2 <- data.frame(standLRT = Exam$standLRT, resid = resid(M2))

qplot(x = standLRT, y = resid, data = true.resid2, geom = c("point", "smooth")) %+% lineup(true = true.resid2, samples = sim.resid2) + facet_wrap( ~ .sample, ncol=4)


# Including standLRT^2 and standLRT^3
M3 <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + (standLRT | school), data = Exam)

qplot(x = standLRT^2, y = resid(M3), data = Exam)

m3.sims  <- simulate(M3, nsim = 19, seed = 1234)
m3.refit <- lapply(m3.sims, refit, object = M3)
m3.resid <- lapply(m3.refit, resid)

sim.resid3 <- do.call("cbind", m3.resid)
sim.resid3 <- melt(sim.resid3)[,-1]
names(sim.resid3) <- c(".n", "resid")
sim.resid3$.n <- as.numeric(str_extract(sim.resid3$.n, "\\d+"))
sim.resid3$standLRT2 <- rep(Exam$standLRT^2, times=19)

true.resid3 <- data.frame(standLRT2 = Exam$standLRT^2, resid = resid(M3))

qplot(x = standLRT2, y = resid, data = true.resid2, geom = "point", alpha = 0.3) %+% lineup(true = true.resid3, samples = sim.resid3) + facet_wrap( ~ .sample, ncol=4)

# Looking at within-school variability
qplot(x = reorder(school, resid(M3), IQR), y = resid(M3), data = Exam, geom = "boxplot") + coord_flip()

sims <- simulate(M3, nsim = 19)

## Refit
refit_fm <- apply(sims, 2, refit, object = M3)

## Extract level-1 residuals
sim_resids <- llply(refit_fm, resid)
sim_resids <- llply(sim_resids, function(x) data.frame(school = Exam$school, resid = x))
true_resids <- data.frame(school = Exam$school, resid = resid(M3))

pushViewport(viewport(layout = grid.layout(5,4)))
vplayout <- function(x, y){viewport(layout.pos.row = x, layout.pos.col = y)}
realp <- qplot(x = reorder(school, resid, IQR), y = resid, data = true_resids, 
               geom = "boxplot", xlab = "school", ylab = "residuals") + coord_flip() + 
                 ylim(-3, 3) + 
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
    p <- qplot(x = reorder(school, resid, IQR), y = resid, 
               data = sim_resids[[j]], geom = "boxplot", xlab = "subject", ylab = "residuals") + 
                 coord_flip() + ylim(-3, 3) + 
                 theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank())
    
    if(!i %in% c(1,5,9,13,17)) p <- p + xlab(NULL)
    if(!i %in% 17:20) p <- p + ylab(NULL)
    
    print(p, vp = vplayout(which(pos == i, arr.ind = TRUE)[1], which(pos == i, arr.ind = TRUE)[2]))
  }
}


# Looking at school-level variables
## Construct school-level data set
SchoolExam <- ddply(Exam, .(school), summarise, size = length(school),
                    schgend = unique(schgend),schavg = unique(schavg),
                    type = unique(type), schLRT = mean(standLRT))

qplot(x = reorder(SchoolExam$schgend, ranef(M3)[[1]][,1], median), 
      y = ranef(M3)[[1]][,1], geom = 'boxplot', 
      xlab = 'school gender', ylab = 'level-2 residual (Intercept)')

qplot(x = reorder(SchoolExam$schgend, ranef(M3)[[1]][,2], median), 
      y = ranef(M3)[[1]][,2], geom = 'boxplot', 
      xlab = 'school gender', ylab = 'level-2 residual (Intercept)')


qplot(x = schavg, y = ranef(M3)[[1]][,1] , data = SchoolExam, 
      geom = c("point", "smooth"), 
      xlab = "average intake score", ylab = "level-2 residual (Intercept)")

qplot(x = schavg, y = ranef(M3)[[1]][,2] , data = SchoolExam, 
      geom = c("point", "smooth"), 
      xlab = "average intake score", ylab = "level-2 residual (Intercept)")


m3.sims  <- simulate(M3, nsim = 19, seed = 1234)
m3.refit <- lapply(m3.sims, refit, object = M3)
m3.resid <- lapply(m3.refit, function(x) ranef(x)[[1]])

sim.resid3 <- do.call("rbind", m3.resid)
sim.resid3$sim <- rownames(sim.resid3)
sim.resid3$.n <- as.numeric(str_extract(sim.resid3$sim, "\\d+"))
sim.resid3$schgend <- rep(SchoolExam$schgend, times=19)
sim.resid3$schavg <- rep(SchoolExam$schavg, times = 19)

true.resid3 <- data.frame(schgend =SchoolExam$schgend, ranef(M3)[[1]])
names(true.resid3) <- c("schgend", "(Intercept)", "standLRT")
true.resid3$schavg <- SchoolExam$schavg

qplot(x = schgend, y =`(Intercept)`, data = true.resid3, geom = "boxplot") %+% lineup(true = true.resid3, samples = sim.resid3) + facet_wrap( ~ .sample, ncol=5)

qplot(x = schavg, y =`(Intercept)`, data = true.resid3, geom = c("point", "smooth"), method = "lm", se=F) %+% lineup(true = true.resid3, samples = sim.resid3) + facet_wrap( ~ .sample, ncol=5)

### Checking for resolution
M4 <- lmer(normexam ~ standLRT + I(standLRT^2) + I(standLRT^3) + sex + schgend + schavg + (standLRT | school), data = Exam)

m4.sims  <- simulate(M4, nsim = 19, seed = 1234)
m4.refit <- lapply(m4.sims, refit, object = M4)
m4.resid <- lapply(m4.refit, function(x) ranef(x)[[1]])

sim.resid4 <- do.call("rbind", m4.resid)
sim.resid4$sim <- rownames(sim.resid4)
sim.resid4$.n <- as.numeric(str_extract(sim.resid4$sim, "\\d+"))
sim.resid4$schgend <- rep(SchoolExam$schgend, times=19)
sim.resid4$schavg <- rep(SchoolExam$schavg, times = 19)

true.resid4 <- data.frame(schgend =SchoolExam$schgend, ranef(M4)[[1]])
names(true.resid4) <- c("schgend", "(Intercept)", "standLRT")
true.resid4$schavg <- SchoolExam$schavg

qplot(x = schgend, y =`(Intercept)`, data = true.resid4, geom = "boxplot") %+% lineup(true = true.resid4, samples = sim.resid4) + facet_wrap( ~ .sample, ncol=5)

qplot(x = schavg, y =`(Intercept)`, data = true.resid4, geom = c("point", "smooth"), method = "lm", se=F) %+% lineup(true = true.resid4, samples = sim.resid4) + facet_wrap( ~ .sample, ncol=5)


### 
M5 <- lmer(normexam ~ standLRT + sex + schgend + (standLRT | school), data = Exam)
qplot(x = standLRT, y = fitted(M5), data = Exam, group = school, geom = "smooth", method = "lm", se=F)

# sims
m5.sims  <- simulate(M5, nsim = 19, seed = 1235)
m5.refit <- lapply(m5.sims, refit, object = M5)
m5.fitted <- lapply(m5.refit, fitted)

sim.fitted5 <- do.call("cbind", m5.fitted)
sim.fitted5 <- melt(sim.fitted5)[,-1]
names(sim.fitted5) <- c(".n", "fitted")
sim.fitted5$.n <- as.numeric(str_extract(sim.fitted5 $.n, "\\d+"))
sim.fitted5$standLRT <- rep(Exam$standLRT, times=19)
sim.fitted5$school <- rep(Exam$school, times=19)

true.fitted5 <- data.frame(fitted = fitted(M5), standLRT = Exam$standLRT, school = Exam$school)

qplot(x = standLRT, y = fitted, data = true.fitted5, group = school, geom = "smooth", method = "lm", se=F) %+% lineup(true = true.fitted5, samples = sim.fitted5) + facet_wrap( ~ .sample, ncol=5)


#-------------------------------------------------------------------------------
# An example of the lines on scatterplots of the random effects that are
# discussed in Morrell and Brant (2000).
#
# Data description: (longituginal data would be best here)
# 
#-------------------------------------------------------------------------------
