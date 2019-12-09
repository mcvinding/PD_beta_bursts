# Logistic regression analysis to get ROC
library(nnet)
library(ROCR)
library(lme4)
library(arm)
library(R.matlab)

# wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

################################################################################
# Functions
# Find optimal cutoff
cut.val <- function(x, mod){
  (logit(x)-coef(mod)[1])/coef(mod)[2]
  # val <- log(x/(1-x)-coef(mod)[1])/coef(mod)[2]
}

################################################################################
## N event analysis
load(file = 'neve.RData')
n.dat1 <- subset(neve.data, session=="1")
n.dat2 <- subset(neve.data, session=="2")

# SESSION 1
# Make logistic model
n1.lmod <- glm(group~nevent, data=n.dat1, family = 'binomial')

# Preditions
pred <- predict(n1.lmod, n.dat1, type="response")
pred <- prediction(pred, n.dat1$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
print(c(Accuaracy=acc, Cutoff=cut))

# ROC
n1.roc <- performance(pred, "tpr", "fpr")

# AUC
n1.auc <- performance(pred, "auc")
n1.auc <- unlist(slot(n1.auc, "y.values"))
print(c(AUC=n1.auc))

# Read cutoff value (not sure this is correct)
n1.cn <- cut.val(cut, n1.lmod)
print(cat("Cutoff: ", round(n1.cn)/3))

# SESSION 23
# Make model
n2.lmod <- glm(group~nevent, data=n.dat2, family = 'binomial')

pred <- predict(n2.lmod, n.dat2, type="response")
pred <- prediction(pred, n.dat2$group)
eval <- performance(pred, "acc")
plot(eval)
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
print(c(Accuaracy=acc, Cutoff=cut))

# ROC
n2.roc <- performance(pred, "tpr", "fpr")

# AUC
n2.auc <- performance(pred, "auc")
n2.auc <- unlist(slot(n2.auc, "y.values"))
print(c(AUC=n2.auc))

# Read cutoff value
n2.cn <- cut.val(cut, n2.lmod)
print(cat("Cutoff: ", round(n2.cn)/3))

################################################################################
## Event duration analysis
load(file = 'leneve.RData')

len_mean <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), mean)
len_medi <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), median)
len_mode <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), Mode)

l.mean.dat1 <- subset(len_mean, session=="1")
l.mean.dat2 <- subset(len_mean, session=="2")

# SESSION 1
# Make logistic model
l1.lmod <- glm(group~x, data=l.mean.dat1, family = 'binomial')

# Preditions
pred <- predict(l1.lmod, l.mean.dat1, type="response")
pred <- prediction(pred, l.mean.dat1$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, l1.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
l1.roc <- performance(pred, "tpr", "fpr")
plot(l1.roc)

# AUC
l1.auc <- performance(pred, "auc")
l1.auc <- unlist(slot(l1.auc, "y.values"))
print(c(AUC=l1.auc))

# SESSION 2
# Make logistic model
l2.lmod <- glm(group~x, data=l.mean.dat2, family = 'binomial')

# Preditions
pred <- predict(l2.lmod, l.mean.dat2, type="response")
pred <- prediction(pred, l.mean.dat2$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, l2.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

l2.cn <- cut.val(cut, l2.lmod)
print(paste("Cutoff: ", l2.cn))

# ROC
l2.roc <- performance(pred, "tpr", "fpr")
plot(l2.roc)

# AUC
l2.auc <- performance(pred, "auc")
l2.auc <- unlist(slot(l2.auc, "y.values"))
print(c(AUC=l2.auc))


################################################################################
## Inter-event interval analysis
load(file = 'itieve.RData')

iti_mean <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), mean)
iti_medi <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), median)
iti_mode <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), Mode)

i.mean.dat1 <- subset(iti_mean, session=="1")
i.mean.dat2 <- subset(iti_mean, session=="2")

# SESSION 1
# Make logistic model
i1.lmod <- glm(group~x, data=i.mean.dat1, family = 'binomial')

# Preditions
pred <- predict(i1.lmod, i.mean.dat1, type="response")
pred <- prediction(pred, i.mean.dat1$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, i1.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
i1.roc <- performance(pred, "tpr", "fpr")
plot(i1.roc)

# AUC
i1.auc <- performance(pred, "auc")
i1.auc <- unlist(slot(i1.auc, "y.values"))
print(c(AUC=i1.auc))

# SESSION 2
# Make logistic model
i2.lmod <- glm(group~x, data=i.mean.dat2, family = 'binomial')

# Preditions
pred <- predict(i2.lmod, i.mean.dat2, type="response")
pred <- prediction(pred, i.mean.dat2$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, i2.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
i2.roc <- performance(pred, "tpr", "fpr")
plot(i2.roc)

# AUC
i2.auc <- performance(pred, "auc")
i2.auc <- unlist(slot(i2.auc, "y.values"))
print(c(AUC=i2.auc))

################################################################################
## Peak amplitude analysis
load(file = 'maxeve.RData')

amp_mean <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), mean)
amp_medi <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), median)
amp_mode <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), Mode)

a.mean.dat1 <- subset(amp_mean, session=="1")
a.mean.dat2 <- subset(amp_mean, session=="2")

# SESSION 1
# Make logistic model
a1.lmod <- glm(group~x, data=a.mean.dat1, family = 'binomial')

# Preditions
pred <- predict(a1.lmod, a.mean.dat1, type="response")
pred <- prediction(pred, a.mean.dat1$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, a1.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
a1.roc <- performance(pred, "tpr", "fpr")
plot(a1.roc)

# AUC
a1.auc <- performance(pred, "auc")
a1.auc <- unlist(slot(a1.auc, "y.values"))
print(c(AUC=a1.auc))

# SESSION 2
# Make logistic model
a2.lmod <- glm(group~x, data=a.mean.dat2, family = 'binomial')

# Preditions
pred <- predict(a2.lmod, a.mean.dat2, type="response")
pred <- prediction(pred, a.mean.dat2$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, a2.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
a2.roc <- performance(pred, "tpr", "fpr")
plot(a2.roc)

# AUC
a2.auc <- performance(pred, "auc")
a2.auc <- unlist(slot(a2.auc, "y.values"))
print(c(AUC=a2.auc))

################################################################################
## Realtive beta power analysis
load(file = 'reldat.RData')

rel.dat1 <- subset(rel.dat, session=="1")
rel.dat2 <- subset(rel.dat, session=="2")

# SESSION 1
# Make logistic model
r1.lmod <- glm(group~relpow, data=rel.dat1, family = 'binomial')

# Preditions
pred <- predict(r1.lmod, rel.dat1, type="response")
pred <- prediction(pred, rel.dat1$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, r1.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
r1.roc <- performance(pred, "tpr", "fpr")
plot(r1.roc)

# AUC
r1.auc <- performance(pred, "auc")
r1.auc <- unlist(slot(r1.auc, "y.values"))
print(c(AUC=r1.auc))

# SESSION 2
#Make logistic model
r2.lmod <- glm(group~relpow, data=rel.dat2, family = 'binomial')

# Preditions
pred <- predict(r2.lmod, rel.dat2, type="response")
pred <- prediction(pred, rel.dat2$group)
eval <- performance(pred, "acc")
plot(eval)                                         # Plot for inspection
max <- which.max(slot(eval, "y.values")[[1]])
acc <- slot(eval, "y.values")[[1]][max]
cut <- slot(eval, "x.values")[[1]][max]
cn <- cut.val(cut, r2.lmod)
print(c(Accuaracy=acc, Cutoff=cut, cut.val=cn))

# ROC
r2.roc <- performance(pred, "tpr", "fpr")
plot(r2.roc)

# AUC
r2.auc <- performance(pred, "auc")
r2.auc <- unlist(slot(r2.auc, "y.values"))
print(c(AUC=r2.auc))

################################################################################

#END