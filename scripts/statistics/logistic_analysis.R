# Logistic regression analysis to get ROC
# library(nnet)
library(ROCR)
library(lme4)
library(arm)
# library(R.matlab)
source("C:/Users/Mikkel/Documents/betabursts/scripts/functions/misc_funs.R")

# wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

################################################################################
# Functions
# Find optimal cutoff
cut.val <- function(x, mod){
  val <- (logit(x)-coef(mod)[1])/coef(mod)[2]
  # val <- (log(x/(1-x))-coef(mod)[1])/coef(mod)[2]
  return(val)
}

# Funtion that do the analysis (implemetet along the way. Clean script and use for all analysis below)
roc_fun <- function(lmod, dat){
  
  # Preditions
  pred <- predict(lmod, dat, type="response")
  pred <- prediction(pred, dat$group)
  eval <- performance(pred, "acc")
  plot(eval)                                         # Plot for inspection
  max <- which.max(slot(eval, "y.values")[[1]])
  acc <- slot(eval, "y.values")[[1]][max]
  cut <- slot(eval, "x.values")[[1]][max]
  cn <- cut.val(cut, lmod)
  print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))
  
  # ROC
  roc <- performance(pred, "tpr", "fpr")
  plot(roc)
  
  # AUC
  auc <- performance(pred, "auc")
  auc <- unlist(slot(auc, "y.values"))
  print(c(AUC=auc))
  
  output <- list("roc" = roc, "auc" = auc, "mod" = lmod)
  return(output)
}

################################################################################
## N event analysis
load(file = 'neve.RData')
n.dat1 <- subset(neve.data, session=="1")
n.dat2 <- subset(neve.data, session=="2")

# SESSION 1
n1.lmod <- glm(group~nevent, data=n.dat1, family = 'binomial')
n1.roc <- roc_fun(n1.lmod, n.dat1)

# SESSION 2
n2.lmod <- glm(group~nevent, data=n.dat2, family = 'binomial')
n2.roc <- roc_fun(n2.lmod, n.dat2)

################################################################################
## Event duration analysis
load(file = 'leneve.RData')

len_mean <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), mean)
len_medi <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), median)
len_mode <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), Mode)

##### MEAN #####
l.mean.dat1 <- subset(len_mean, session=="1")
l.mean.dat2 <- subset(len_mean, session=="2")

# SESSION 1
l1mean.lmod <- glm(group~x, data=l.mean.dat1, family = 'binomial')
l1mean.roc <- roc_fun(l1mean.lmod, l.mean.dat1)

# SESSION 2
l2mean.lmod <- glm(group~x, data=l.mean.dat2, family = 'binomial')
l2mean.roc <- roc_fun(l2mean.lmod, l.mean.dat2)

##### MEDIAN #####
l.medi.dat1 <- subset(len_medi, session=="1")
l.medi.dat2 <- subset(len_medi, session=="2")

# SESSION 1
l1medi.lmod <- glm(group~x, data=l.medi.dat1, family = 'binomial')
l1medi.roc <- roc_fun(l1mean.lmod, l.medi.dat1)

# SESSION 2
l2medi.lmod <- glm(group~x, data=l.medi.dat2, family = 'binomial')
l2medi.roc <- roc_fun(l2medi.lmod, l.medi.dat2)

##### MODE #####
l.mode.dat1 <- subset(len_mode, session=="1")
l.mode.dat2 <- subset(len_mode, session=="2")

# SESSION 1
l1mode.lmod <- glm(group~x, data=l.mode.dat1, family = 'binomial')
l1mode.roc <- roc_fun(l1mode.lmod, l.mode.dat1)

# SESSION 2
l2mode.lmod <- glm(group~x, data=l.mode.dat2, family = 'binomial')
l2mode.roc <- roc_fun(l2mode.lmod, l.mode.dat2)

################################################################################
## Inter-event interval analysis
load(file = 'itieve.RData')

iti_mean <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), mean)
iti_medi <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), median)
iti_mode <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), Mode)

##### MEAN #####
i.mean.dat1 <- subset(iti_mean, session=="1")
i.mean.dat2 <- subset(iti_mean, session=="2")

# SESSION 1
i1mean.lmod <- glm(group~x, data=i.mean.dat1, family = 'binomial')
i1mean.roc <- roc_fun(i1mean.lmod, i.mean.dat1)

# SESSION 2
i2mean.lmod <- glm(group~x, data=i.mean.dat2, family = 'binomial')
i2mean.roc <- roc_fun(i2mean.lmod, i.mean.dat2)

##### MEDIAN #####
i.medi.dat1 <- subset(iti_medi, session=="1")
i.medi.dat2 <- subset(iti_medi, session=="2")

# SESSION 1
i1medi.lmod <- glm(group~x, data=i.medi.dat1, family = 'binomial')
i1medi.roc <- roc_fun(i1medi.lmod, i.medi.dat1)

# SESSION 2
i2medi.lmod <- glm(group~x, data=i.medi.dat2, family = 'binomial')
i2medi.roc <- roc_fun(i2medi.lmod, i.medi.dat2)

##### MODE #####
i.mode.dat1 <- subset(iti_mode, session=="1")
i.mode.dat2 <- subset(iti_mode, session=="2")

# SESSION 1
# Make logistic model
i1mode.lmod <- glm(group~x, data=i.mode.dat1, family = 'binomial')
i1mode.roc <- roc_fun(i1mode.lmod, i.mode.dat1)

# SESSION 2
i2mode.lmod <- glm(group~x, data=i.mode.dat2, family = 'binomial')
i2mode.roc <- roc_fun(i2mode.lmod, i.mode.dat2)

################################################################################
## Peak amplitude analysis
load(file = 'maxeve.RData')

amp_mean <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), mean)
amp_medi <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), median)
amp_mode <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), Mode, digits=2)

a.mean.dat1 <- subset(amp_mean, session=="1")
a.mean.dat2 <- subset(amp_mean, session=="2")

# session 1
a1mean.lmod <- glm(group~x, data=a.mean.dat1, family = 'binomial')
a1mean.roc <- roc_fun(a1mean.lmod,a.mean.dat1)

# SESSION 2
a2mean.lmod <- glm(group~x, data=a.mean.dat2, family = 'binomial')
a2mean.roc <- roc_fun(a2mean.lmod,a.mean.dat2)

##### MEDIAN #####
a.medi.dat1 <- subset(amp_medi, session=="1")
a.medi.dat2 <- subset(amp_medi, session=="2")

# session 1
a1medi.lmod <- glm(group~x, data=a.medi.dat1, family = 'binomial')
a1medi.roc <- roc_fun(a1medi.lmod,a.medi.dat1)

# SESSION 2
a2medi.lmod <- glm(group~x, data=a.medi.dat2, family = 'binomial')
a2medi.roc <- roc_fun(a2medi.lmod,a.medi.dat2)

##### MODE #####
a.mode.dat1 <- subset(amp_mode, session=="1")
a.mode.dat2 <- subset(amp_mode, session=="2")

# session 1
a1mode.lmod <- glm(group~x, data=a.mode.dat1, family = 'binomial')
a1mode.roc <- roc_fun(a1mode.lmod,a.mode.dat1)

# SESSION 2
a2mode.lmod <- glm(group~x, data=a.mode.dat2, family = 'binomial')
a2mode.roc <- roc_fun(a2mode.lmod,a.mode.dat2)

################################################################################
## Realtive beta power analysis
load(file = 'reldat.RData')

rel.dat1 <- subset(rel.dat, session=="1")
rel.dat2 <- subset(rel.dat, session=="2")

# SESSION 1
r1.lmod <- glm(group~relpow, data=rel.dat1, family = 'binomial')
r1.roc <- roc_fun(r1.lmod, rel.dat1)

# SESSION 2
r2.lmod <- glm(group~relpow, data=rel.dat2, family = 'binomial')
r2.roc <- roc_fun(r2.lmod, rel.dat2)

################################################################################
# 

################################################################################
# Save
save("r1.roc", "r2.roc", 
     "n1.roc", "n2.roc",
     "l1mean.roc", "l2mean.roc", "l1medi.roc", "l2medi.roc", "l1mode.roc", "l2mode.roc",
     "i1mean.roc", "i2mean.roc", "i1medi.roc", "i2medi.roc", "i1mode.roc", "i2mode.roc",
     "a1mean.roc", "a2mean.roc", "a1medi.roc", "a2medi.roc", "a1mode.roc", "a2mode.roc",
     file="C:/Users/Mikkel/Documents/betabursts/groupanalysis/roc_results.Rdata")

#END