source("R/add_interaction.R")
library(grid)
library(ggplot2)


##### lineup number 1 #####
load("exam-fanned.RData")

location <- 16
ggplot(aes(x = standLRT, y = y, group=school), data=gcse) + geom_smooth(method="lm", se=F, colour=rgb(0, 0, 0, alpha=0.5),  alpha=0.1) + facet_wrap(~sample) +
	theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
ggsave("normexam_fanned_lineup13.pdf")

make_interactive(filename= sprintf("exam-fanned-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("exam-fanned-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")

##### lineup number 2 #####

load("autism-fanned2.RData")

location <- 1
qplot(x = age.2, y = y, data = autism.true.y, group = childid, geom = "line", se=F, alpha = I(0.3)) %+% lineup(true = autism.true.y, samples = mod2.sim.y, pos=location) + facet_wrap( ~ .sample, ncol=5) +
theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
                 axis.text.x = element_blank(), panel.grid.major.y = element_blank(), 
                 axis.title=element_blank()) 
ggsave(sprintf("figures/autism-fanned2-%s-multiple.svg", location)

make_interactive(filename= sprintf("autism-fanned2-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("autism-fanned2-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="toggle")




##### lineup number 3 #####

location <- 1
qplot(x = sicdegp, y = residual, data = M2.true.sicdegp, geom = "boxplot", 
	fill = sicdegp, outlier.size = 2, alpha=I(0.6)) %+% 
	lineup(true = M2.true.sicdegp, samples = M2.sim.sicdegp, pos=location) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10) + 
	ylab("level-1 residuals") + 
	scale_fill_brewer("", palette="Set2", labels=c("low", "medium", "high")) +
	theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
ggsave("autism_sicdegp_level1_lineup5.pdf")
save(M2.true.sicdegp, M2.sim.sicdegp, file="autism-ordered.RData")

make_interactive(filename= sprintf("autism-ordered-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("autism-ordered-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")


##### lineup number 4 #####

location <- 9

qplot(x = sicdegp, y = residual, data = M2.true.sicdegp, geom = "boxplot", 
	fill = sicdegp, outlier.size = 2, alpha=I(0.6)) %+% 
	lineup(true = M2.true.sicdegp, samples = M2.sim.sicdegp) + 
	facet_wrap( ~ .sample, ncol=5) + 
	ylim(-10, 10) + 
	ylab("level-1 residuals") + 
	scale_fill_brewer("", palette="Set2", labels=c("high", "low", "medium")) +
	theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
ggsave("autism_sicdegp_level1_lineup4.pdf")
save(M2.true.sicdegp, M2.sim.sicdegp, file="autism-unordered.RData")

make_interactive(filename= sprintf("autism-unordered-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("autism-unordered-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")

##### lineup number 5 #####
load("exam-corr.RData")
location <- 7

qplot(x = `(Intercept)`, y = standLRT, data = true.M2.ranef, 
	geom = c("point", "smooth"), method = "lm", se = F, alpha = I(0.4)) %+% 
	lineup(true = true.M2.ranef, samples = M3.sim.ranef) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
ggsave("exam_corr_lineup7.pdf")

make_interactive(filename= sprintf("exam-corr-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("exam-corr-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")


##### lineup number 6 #####

load("dialyzer-nonlinear.RData")

location <- 17
qplot(pressure, resid, data = m1.resid.df, geom = c("point", "smooth"),
	method = "loess") %+%
  	lineup(true = m1.resid.df, samples = m1.sim.resids, pos=location) +
  	facet_wrap(~ .sample, ncol = 5) +
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

make_interactive(filename= sprintf("dialyzer-nonlinear-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("dialyzer-nonlinear-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")




##### lineup number 7 #####
ggsave("exam-homogeneity.RData")
location <- 7

qplot(standLRT2, resid, data = lev1_resid_fm2,
      	geom = "point", alpha = I(0.4)) %+%
  	lineup(true = lev1_resid_fm2, samples = sim_fm2_lev1_resid) +
  	facet_wrap(~ .sample, ncol = 5) +
  	geom_hline(aes(yintercept = 0), colour = I("white"), alpha = I(0.5)) + 
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
make_interactive(filename= sprintf("exam-homogeneity-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("exam-homogeneity-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")




##### lineup number 8 #####
load(file="cyclone-bad.RData")

qplot(x = factor(rank), y = EB.resid, data = resids, 
               geom = "boxplot", xlab = "", ylab = "", outlier.size = 1.5) + coord_flip() + 
                 ylim(-150, 150) + 
                 theme(plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.y = element_blank(),
                 axis.text.x = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank()) +
				facet_wrap(~sample)

location <- resids$sample[nrow(resids)]
make_interactive(filename= sprintf("cyclone-bad-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("cyclone-bad-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="toggle")

##### lineup number 9 #####
load("radon-tsim-6.RData")
location <- 6
# Lineup of random slopes
qplot(sample = basement, data = b1, stat = "qq") %+%
	lineup(true = b1, sample = sim.b1, pos=location) + 
	facet_wrap(~ .sample, ncol = 5) + 
	geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
#	xlab("Normal Quantiles") + ylab("Sample Quantiles") +
	ylab(NULL) + xlab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

make_interactive(filename= sprintf("radon-tsim-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("radon-tsim-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")

##### lineup number 10 #####
load("dialyzer-nonlinear.RData")


location <- 3
qplot(pressure, resid, data = m2.resid.df,
      	geom = "point", alpha = I(0.5)) %+%
  	lineup(true = m2.resid.df, samples = m2.sim.resids) +
  	facet_wrap(~ .sample, ncol = 5) +
 	 xlab(NULL) + ylab(NULL) + 
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
	axis.ticks.x = element_blank(), axis.ticks.y = element_blank())

make_interactive(filename= sprintf("dialyzer-heterogeneous-%s-multiple.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js")
make_interactive(filename= sprintf("dialyzer-heterogeneous-%s-single.svg", location), 
		script="http://www.hofroe.net/examples/lineup/action.js", toggle="select")
