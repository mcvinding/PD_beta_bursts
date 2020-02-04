# Between group-session analysis of relative beta power
library(R.matlab)
library(BayesFactor)
library(brms)
library(plyr)

# Define paths
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

# Load data
temp <- readMat("B_relpow.mat")
ptns.relpow1 <- temp$ptns.relpow1
ptns.relpow2 <- temp$ptns.relpow2
ctrl.relpow1 <- temp$ctrl.relpow1
ctrl.relpow2 <- temp$ctrl.relpow2

# Arrange data
relpow <- c(ptns.relpow1, ptns.relpow2, ctrl.relpow1, ctrl.relpow2)
subs <- as.factor(c(rep(1:length(ptns.relpow1),2), rep(1:length(ctrl.relpow1)+20,2)))
session <- as.factor(c(rep(c(1,2), each=length(ptns.relpow1)), rep(c(1,2), each=length(ctrl.relpow1)*2)))
group <- as.factor(c(rep(1, length(ptns.relpow1)*2), rep(2, length(ctrl.relpow1)*2)))

rel.dat <- data.frame(relpow = c(ptns.relpow1, ptns.relpow2, ctrl.relpow1, ctrl.relpow2),
                      subj = as.factor(c(rep(1:length(ptns.relpow1),2), rep(1:length(ctrl.relpow1)+20,2))),
                      session = as.factor(c(rep(c(1,2), each=length(ptns.relpow1)), rep(c(1,2), each=length(ctrl.relpow1)))),
                      group = as.factor(c(rep(1, length(ptns.relpow1)*2), rep(2, length(ctrl.relpow1)*2)))
                      )

rel.dat$group <- revalue(rel.dat$group, c("1"="ptns", "2"="ctrl"))

# Save data
save(rel.dat, file = 'reldat.RData')

###################################################################
# Summary
b.summary <- aggregate(rel.dat$relpow, list(rel.dat$group, rel.dat$session), mean)
names(b.summary) <- c("Group","Session","mean")
b.summary.sd <- aggregate(rel.dat$relpow, list(rel.dat$group, rel.dat$session), sd)
b.summary$sd <- b.summary.sd$x

###################################################################
# ANALYSIS
# Bayesian "T-tests"
ttestBF(ptns.relpow1, ptns.relpow2, paird = TRUE)
ttestBF(ctrl.relpow1, ctrl.relpow2, paird = TRUE)
ttestBF(ptns.relpow1, ctrl.relpow1, paird = FALSE)
ttestBF(ptns.relpow2, ctrl.relpow2, paird = FALSE)

ptns.diff <- ptns.relpow1-ptns.relpow2
ctrl.diff <- ctrl.relpow1-ctrl.relpow2

ttestBF(ptns.diff, ctrl.diff, paird = FALSE)

# Student's T-test
t.test(ptns.relpow1, ptns.relpow2, paird = TRUE)
t.test(ctrl.relpow1, ctrl.relpow2, paird = TRUE)
t.test(ptns.relpow1, ctrl.relpow1, paird = FALSE)
t.test(ptns.relpow2, ctrl.relpow2, paird = FALSE)

t.test(ptns.diff, ctrl.diff, paird = FALSE)
