setwd("~/papers/2014-sociology_chapter/")
library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(mlmRev)   # for the exam data
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(gridSVG)

make_interactive <- function(filename, script, toggle="toggle", background="#e5e5e5") {
# use #ffffff for white background
  require(gridSVG)
  grobs <- grid.ls()
  
  idx <- grep("panel-", grobs$name)
  for (i in idx) { 
    grid.garnish(grobs$name[i],
                 onmouseover=paste("frame('",grobs$name[i+2], ".1')", sep=""),
                 onmouseout=paste("deframe('",grobs$name[i+2], ".1')", sep=""), 
                 onmousedown=sprintf("%shigh(evt, '%s.1', '%s')", toggle, grobs$name[i+2], background)
                 )
  }

  # use script on server to get locally executable javascript code
  # or use inline option
  grid.script(filename=script)
  grid.export(filename, uniqueNames=FALSE, exportJS="inline", exportCoords="inline", exportMappings="inline")
}

#################### exam-fanned-withslope #####################
###### create data

set.seed("20140501")
for (rep in 1:5) {
## Adding the random slope
M2 <- lmer(normexam ~ standLRT + (standLRT  | school), data = Exam)

M2.sims  <- simulate(M2, nsim = 19)#, seed = 1234)
M2.refit <- lapply(M2.sims, refit, object = M2)
M2.simy <- lapply(M2.refit, function(x) x@resp$y)

sim2.y <- do.call("cbind", M2.simy)
sim2.y <- melt(sim2.y)[,-1]
names(sim2.y) <- c(".n", "y")
sim2.y$.n <- as.numeric(str_extract(sim2.y$.n, "\\d+"))
sim2.y$standLRT <- rep(Exam$standLRT, rep = 19)
sim2.y$school <- rep(Exam$school, rep = 19)

# The "harder way" to get alpha working
M2.fitted <- ddply(sim2.y, .(.n, school), function(x) {
	m <- lm(y ~ standLRT, data = x)
	data.frame(x, fitted = fitted(m))
})

true.df <- subset(Exam, select = c(school, normexam, standLRT))
colnames(true.df)[2] <- "y"

true.fitted <- ddply(true.df, .(school), function(x) {
	m <- lm(y ~ standLRT, data = x)
	data.frame(x, fitted = fitted(m))	
})


# save(M2.fitted, true.fitted, file = "exam-fanned-withslope.RData")

datfile <- rbind(M2.fitted, data.frame(.n=20, true.fitted))
resample <- sample(20)
location <- resample[20]
datfile$.sample <- resample[datfile$.n]
datfile$true <- location

file_name <- sprintf("~/turk14/lineups/data/exam-fanned-withslope-%d-%d.csv",rep,location)
write.csv(datfile[,-1],file=file_name, row.names=FALSE)
}

#################### exam-fanned-withslope #####################
###### create plot
setwd("~/turk14/lineups/data/")
datadir <- dir()

files <- datadir[grep("exam-fanned-withslope",datadir)]
for (fname in files) {
	dframe <- read.csv(fname)
	qplot(x = standLRT, y = fitted, data = dframe, group = school, 
		geom = "line", alpha = I(0.5)) +
		facet_wrap( ~ .sample, ncol=5) + 
		xlab(NULL) + ylab(NULL) +
		theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
			axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

	ggsave(file=gsub("csv","pdf", fname))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[4]
	location <- params[5]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("exam.fanned.with.slope-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=fname,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)

	}


}

#####################################################

#################### dialyzer-heterogeneous #####################
# data copied from turk 11
library(RColorBrewer)

###### create plot
setwd("~/turk14/lineups/data/")
datadir <- dir()

files <- datadir[grep("dialyzerheterogeneous",datadir)]
for (fname in files) {
	dframe <- read.csv(fname)
	dframe$pressurecat <- cut(dframe$pressure,breaks=c(0,0.45,0.75,1.25,1.75, 2.25, 2.75, 3.25))

	print(qplot(pressurecat, resid, data = dframe, geom = "boxplot", fill=pressurecat, alpha=I(0.8)) +
  	facet_wrap(~ sample, ncol = 5) +
 	 xlab(NULL) + ylab(NULL) + scale_fill_brewer() +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position="none") + 
	scale_x_discrete(expand = c(.1,.1))
	)

	ggsave(file=gsub("csv","pdf", fname))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[2]
	location <- params[3]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("dialyzer.heterogeneous-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=fname,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)

	}


}

###################################
######## dialyzer-homogeneous ######

##################################
# create data


library(MEMSS)
M2 <- lmer(rate ~ (pressure + I(pressure^2) + I(pressure^3) + I(pressure^4))*QB + (pressure + I(pressure^2) | Subject), data = Dialyzer)

m2.resid.df <- data.frame(.n=20, Subject = M2@frame$Subject, resid = resid(M2))


set.seed("20040110")
foo <- rnorm(150000)
rm(foo)

for (null in 1:5) {
m2.sims <- simulate(M2, nsim = 19)
m2.refit <- lapply(m2.sims, refit, object = M2)
m2.sim.resids <- ldply(m2.refit, function(x) data.frame(x@frame, resid = resid(x)))
m2.sim.resids$.n <- as.numeric(str_extract(m2.sim.resids$.id, "\\d+"))

m2.resid <- rbind(m2.sim.resids[,c(".n", "Subject", "resid")], m2.resid.df)
m2.resid$sample <- sample(20,20, replace=FALSE)[m2.resid$.n]
location <- m2.resid$sample[nrow(m2.resid)]
write.csv(m2.resid, sprintf("dialyzerhomogeneous-%s-%s.csv", null, location), row.names=FALSE)
}


##################################
# create plot - boxplots

setwd("~/turk14/lineups/data/")
datadir <- dir()

files <- datadir[grep("dialyzerhomogeneous",datadir)]
for (fname in files) {
	dframe <- read.csv(fname)

	print(qplot(Subject, resid, data = dframe, geom = "boxplot", fill=I("grey70"), alpha=I(0.8)) +
  	facet_wrap(~ sample, ncol = 5) +
#  	geom_hline(yintercept=0, colour="grey50") +
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position="none") + 
	scale_x_discrete(expand = c(.075,.075))
	)

	ggsave(file=gsub("csv","pdf", fname))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[2]
	location <- params[3]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("dialyzer.homogeneous-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=fname,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)

	}


}

##################################
# create plot - dotplots

setwd("~/turk14/lineups/data/")
datadir <- dir()

files <- datadir[grep("dialyzerhomogeneous",datadir)]
for (fname in files) {
	dframe <- read.csv(fname)

dframe <- ddply(dframe, .(sample, Subject), transform, var = var(resid))
dframe <- ddply(dframe, .(sample), transform, subid = rank(var, ties.method="min"))

	print(qplot(subid, resid, data = dframe, geom = c("point"), fill=I("grey70"), alpha=I(0.5)) +
  	facet_wrap(~ sample, ncol = 5) +
  	geom_hline(yintercept=0, colour="red") +
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position="none") + 
	scale_x_continuous(expand = c(.075,.075))
	)
	fname <- gsub("ous", "ous-dots", fname, fixed=TRUE)
	ggsave(file=gsub("csv","pdf", fname))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[2]
	location <- params[3]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("dialyzer.homogeneous.dots-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=fname,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)

	}


}

###################################
######## cyclone plots for radon data ######

##################################
# create data

#-------------------------------------------------------------------------------
# Script for constructing the lineups of cyclone plots for the radon data
#
# Adam Loy, Heike Hofmann, Dianne Cook
# June 2013
#-------------------------------------------------------------------------------

library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(grid)

radon <- ddply(radon, .(county), transform, n = length(county))

# The model - with fewer counties
fm <- lmer(log.radon ~ basement + uranium + (basement | county), data = subset(radon, n > 4))

# The *observed residuals* and other variables needed for lineup
#observed.cyclone.df <- data.frame(fm@frame, resid = resid(fm))
observed.cyclone.df <- HLMresid(fm, type="LS", level=1)
set.seed(987654321)

for (null in 1:5) {
# Simulating a set of 19 null plots
sim.y   <- simulate(fm, nsim = 19) #, seed = 987654321) ## Change seed for diff. nulls                       
sim.mod <- apply(sim.y, 2, refit, object = fm)

#null.cyclone.df <- ldply(sim.mod, function(x) data.frame(x@frame, resid = resid(x)))
null.cyclone.df <- ldply(sim.mod, function(x) HLMresid(x, type="LS", level=1))
null.cyclone.df$.id <- as.numeric( str_extract(null.cyclone.df$.id, "\\d+") )

cyclone.df <- rbind.fill(observed.cyclone.df, null.cyclone.df)

plot.order <- sample.int(20, 20)
location <- plot.order[1]
cyclone.df$.id <- rep(plot.order, each = nrow(fm@frame))

cyclone.df$resid <- cyclone.df$LS.resid
cyclone.df <- ddply(cyclone.df, .(.id, county), transform, iqr = IQR(resid, na.rm=TRUE))
cyclone.df <- ddply(cyclone.df, .(.id), transform, rank = order(order(iqr, county)))
cyclone.df <- ddply(cyclone.df, .(.id, county), transform, rank = min(rank))
cyclone.df <- dlply(cyclone.df, .(.id), transform, rank = rank(rank))
cyclone.df <- ldply(cyclone.df, function(df) {
	df$rank <- factor(df$rank)
	levels(df$rank) <- 1:length(levels(df$rank))
	return(df)
})
cyclone.df$location <- location
fname <- sprintf("data/cyclone.radon-%s-%s.csv", null, location)

write.csv(cyclone.df, file=fname, row.names=FALSE)
}

##################################
# create plot

setwd("~/turk14/lineups/data/")
datadir <- dir()

files <- datadir[grep("cyclone.radon",datadir)]
for (fname in files) {
	dframe <- read.csv(fname)


	print(qplot(x = factor(rank), y = resid, data = dframe, geom = "boxplot", outlier.size = 1.5) + 
	coord_flip() + 
	ylim(-2, 2) + 
	xlab(NULL) + 
	ylab(NULL) + 
	theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
	      axis.text.x = element_blank(), axis.ticks.y = element_blank(),
	      axis.ticks.x = element_blank(), 
	      panel.grid.major.y = element_blank()) +
	facet_wrap(~ .id)
	)

	ggsave(file=gsub("csv","pdf", fname))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[2]
	location <- params[3]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("radon.cyclone-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=fname,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)

	}
}


###################################
######## Q-Q plots for slope random effect in radon data ######

##################################
# prepare model

### Preliminaries
library(lme4)
library(arm)       # for nicer print outs from lmer
library(ggplot2)
library(nortest)   # for tests of normality
library(boot)      # for simulation envelopes

### Reading in the radon data provided by Gelman and Hill (2007)
srrs2 <- read.table ("~/papers/2013-minconfounded_chapter/data/srrs2.dat", header=T, sep=",")

## Restricting attention to Minnesota
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
basement <- floor

## Getting county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

## Reading the county-level data
srrs2.fips <- srrs2$stfips*1000 + srrs2$cntyfips
cty <- read.table ("~/papers/2013-minconfounded_chapter/data/cty.dat", header=T, sep=",")
usa.fips <- 1000*cty[,"stfips"] + cty[,"ctfips"]
usa.rows <- match (unique(srrs2.fips[mn]), usa.fips)
uranium <- cty[usa.rows,"Uppm"]
u <- log (uranium)
uranium <- u[county]

library(HLMdiag)
sim_env <- function(x, conf = .95, line=FALSE){
  n <- length(x)
  P <- ppoints(x)[rank(x)]
  z <- qnorm(P)
  if (line) {
    a <- as.numeric(HLMdiag:::qqlineInfo(x)[1])
    b <- as.numeric(HLMdiag:::qqlineInfo(x)[2])
  } else {
    a <- 0  # we know that the line should be the identity
    b <- 1
  }
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- a + b * z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
  return(data.frame(lower, upper, fit.value))
  
}
### Fitting a model with random intercept for county and 
### random slope for basement.
fm <- lmer(log.radon ~ basement + uranium + (basement | county), REML = FALSE)

#############################
# create data
setwd("~/turk14/lineups/data/")


getneweffects <- function(model) {
	y.sim <- simulate(model, nsim = 1)	# simulate responses
	m <- refit(model, newresp = y.sim)
	b <- ranef(m)[[1]]
	names(b) <- c("b1", "b2")
	if (diff(range(b$b1)) == 0 | diff(range(b$b2)) == 0) return(getneweffects(model))
	else return(m)
}

set.seed(20140501)
for (null in 1:5) {
	mod.sim <- llply(1:19, function(x) getneweffects(model=fm))

	mod.sim[[20]] <- fm
	names(mod.sim) <- 1:20
	
	# get random effects out:
	bs <- ldply(mod.sim, function(x) ranef(x)[[1]])
	names(bs) <- c("id", "b1", "b2")

	sd1 <- HLMdiag:::qqlineInfo(bs$b1[bs$id==20])[2]
	sd2 <- HLMdiag:::qqlineInfo(bs$b2[bs$id==20])[2]
	bs$b1 <- bs$b1/sd1
	bs$b2 <- bs$b2/sd2
	
	ds <- ddply(bs, .(id), transform, 
		b1.qq=qqnorm(b1, plot.it=FALSE), 
		b2.qq=qqnorm(b2, plot.it=FALSE))
	envs <- ddply(ds, .(id), transform, 
		b1.env=sim_env(b1), 
		b2.env=sim_env(b2))	

	samples <- sample(20)
	envs$samples <- samples[as.numeric(envs$id)]
	location <- samples[20]
	envs$plot_location <- location
	fname <- sprintf("radon.qq-%d-%d.csv", null, location)
	
	write.csv(envs, file=fname, row.names=FALSE)
}


#########################
# Q-Q plots

# create plots

files <- dir()
files <- files[grep("radon.qq", files)]

for (fname in files) {
	dframe <- read.csv(fname)

	
	##### random effects b1
	b1name <- gsub("qq-", "qq.b1-", fname)
		
	print(
	ggplot(data=dframe, aes(x=b1.qq.x, y=b1.qq.y)) + 
	  geom_ribbon(aes(x=b1.qq.x, ymin=b1.env.lower, ymax=b1.env.upper), alpha=0.2) +
	  geom_abline(colour = I("grey60")) + 
      geom_point() + 
      facet_wrap(~samples, ncol=5) +  xlab("") + ylab("")+
      theme(axis.text=element_blank(), axis.ticks=element_blank(),    
	      plot.margin = unit(c(.1,.1,.1,.1), "cm")) 
	)

	ggsave(file=gsub("csv","pdf", b1name))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[2]
	location <- params[3]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("radon.qq-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=b1name,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)
	}

	##### random effects b2
	b2name <- gsub("qq-", "qq.b2-", fname)
	
	print(
	ggplot(data=dframe, aes(x=b2.qq.x, y=b2.qq.y)) + 
	  geom_ribbon(aes(x=b2.qq.x, ymin=b2.env.lower, ymax=b2.env.upper), alpha=0.2) +
	  geom_abline(colour = I("grey60")) + 
      geom_point() + 
      facet_wrap(~samples, ncol=5) +  xlab("") + ylab("")+
      theme(axis.text=element_blank(), axis.ticks=element_blank(),    
	      plot.margin = unit(c(.1,.1,.1,.1), "cm")) 
	)


	ggsave(file=gsub("csv","pdf", b2name))
	tmpfile <- tempfile(tmpdir="~")
	params <- unlist(strsplit(gsub(".csv", "", fname),"-"))
	rep <- params[2]
	location <- params[3]
	for (click in c("multiple", "single")) {
		pic_name <- sprintf("%s-%s.svg", gsub("~/","", tmpfile), click)
		sample_size <- nrow(dframe)/20

		test_param <- sprintf("turk14-%s", click)
		param_value <- sprintf("radon.qq-%s", rep)

		write.table(data.frame(
				sample_size=sample_size, 
				test_param=test_param,
				param_value=param_value,
				p_value=NA,
				obs_plot_location=location, 
				pic_name=pic_name,
				experiment="turk14",
				difficulty=b2name,
				data_name=fname
				), 
			file="../picture-details.csv", row.names=FALSE, sep=",",
			col.names=!file.exists("../picture-details.csv"), append=TRUE)

		toggle = "toggle"
		if (click =="single") toggle="select"
		make_interactive(filename= pic_name, 
			script="http://www.hofroe.net/examples/lineup/action.js", toggle=toggle)
	}
}

