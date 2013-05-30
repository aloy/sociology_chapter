source("R/add_interaction.R")

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


##### lineup number 3 #####

location <- 10
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



