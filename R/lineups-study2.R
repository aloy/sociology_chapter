
library(grid)
library(ggplot2)

##### lineup number 1 #####
load("exam-fanned-withslope.RData")

location <- 18
qplot(x = standLRT, y = fitted, data = true.fitted, group = school, 
	geom = "line", alpha = I(0.5)) %+% 
	lineup(true = true.fitted, samples = M2.fitted, pos = 18) + 
	facet_wrap( ~ .sample, ncol=5) + 
	xlab(NULL) + ylab(NULL) +
	theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
		axis.ticks.x = element_blank(), axis.ticks.y = element_blank())


