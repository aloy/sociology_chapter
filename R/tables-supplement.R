require(plyr)
require(magrittr)
require(vinference)
require(reshape2)
require(xtable)

set.seed(20161001)
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
