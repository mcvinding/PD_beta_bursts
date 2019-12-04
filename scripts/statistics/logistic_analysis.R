# Logistic regression analysis to get ROC
library(nnet)
library(ROCR)
library(lme4)

# wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

## N event analysis
################################################################################
load(file = 'neve.RData')
load(file = 'leneve.RData')
load(file = 'itieve.RData')
load(file = 'reldat.RData')


str(neve.data)
str(itieve.data)
str(rel.dat)

# Playing with data
dat <- subset(neve.data, session=="1")

sum.dat <- aggregate(eve.iti ~ group+session+subs, itieve.data, median)
dat <- subset(rel.dat, session=="1")

# Logistic regression Model
lmod2 <- multinom(group~nevent, data=dat)
lmod <- glm(group~nevent, data=dat, family = 'binomial')
lmod <- glm(group~eve.iti, data=dat, family = 'binomial')
lmod <- glm(group~relpow, data=dat, family = 'binomial')

# Table (currently not working for glmer mods)
p <- predict(lmod, dat)
tab <- table(p, dat$group)
sum(diag(tab))/sum(tab)
1-sum(diag(tab))/sum(tab)

# Plot regression (this behaves wierd - look at old plot to see why)
idx <- ifelse(dat$group==2,1,0)
newdat <- data.frame(relpow=seq(min(dat$relpow), max(dat$relpow),len=100))
newdat$p = predict(lmod, newdata=newdat, type="response")

plot(dat$relpow, idx)
plot(newdat$p ~ relpow, newdat, lwd=2)

# Model performance
# pred <- predict(lmod, dat, type="prob")   # with mulitnom
pred <- predict(lmod, dat, type="response")
pred <- prediction(pred, dat$group)
eval <- performance(pred, "acc")
plot(eval)
# abline(h=.84)
# abline(v=.41)

# Best value
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
print(c(Accuaracy=acc, Cutoff=cut))

# ROC
roc <- performance(pred, "tpr", "fpr")
plot(roc,
     colorize=T)
abline(a=0, b=1)

# AUC
auc <- performance(pred, "auc")
unlist(slot(auc, "y.values"))

# Read cutoff value (not sure this is correct)
cn <- (invlogit(cut)-coef(lmod)[1])/coef(lmod)[2]
abline(v=cn)
